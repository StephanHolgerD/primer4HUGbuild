'''
streamlit run primer4/stream.py -- --config config.json
pytest --config config.json
'''


import argparse
import json
from pathlib import Path
import pdb
import os

import click
from cdot.hgvs.dataproviders import JSONDataProvider
import gffutils
import hgvs
import pandas as pd
from pyfaidx import Fasta
import streamlit as st

from primer4.models import Variant, ExonDelta, SingleExon, ExonSpread, Template
from primer4.design import design_primers, check_for_multiple_amplicons
from primer4.hacks import download_button
from primer4.utils import (
    log,
    mask_sequence,
    reconstruct_mrna,
    sync_tx_with_feature_db,
    )
from primer4.vis import prepare_data_for_vis, primers_to_df
# from primer4.vis import Beauty, prepare_mock_data_for_vis
from primer4.warnings import warn


def gimme_some_primers(method, code, fp_genome, genome, hdp, db, vardbs, params, max_variation):
    '''
    gimme_some_primers('sanger', 'NM_000546.6:c.215C>G', ...)
    gimme_some_primers('qpcr', ('NM_001145408.2', 6), ...)
    gimme_some_primers('mrna', ('NM_000546.6', 6, 7), ...)
    '''
    if method == 'sanger':
        try:
            v = Variant(code[0], hdp, db)
        # TODO: Catch this error in the main fn so we don't have st.x() fn
        # all over the code? We could let all errors propagate up to there ...
        except hgvs.exceptions.HGVSParseError:
            st.warning('Variant could not be parsed, check syntax (spaces?).')
            st.stop()

    elif method == 'qpcr':
        v = SingleExon(code[0], int(code[1]))
    
    elif method == 'mrna':
        v = ExonSpread(code[0], int(code[1]), int(code[2]))
    
    else:
        raise ValueError('Method is not implemented, exit.')

    tmp = Template(v, db)
    
    # Annotate templete
    tmp.load_variation_freqs_(vardbs)
    tmp.load_variation_(max_variation)
    tmp.get_sequence_(genome)

    # Mask and get primers
    if method == 'sanger':
        masked = mask_sequence(tmp.sequence, tmp.mask)
        constraints = tmp.apply(method, db, params)
        primers = [p for p in next(design_primers(masked, constraints, params, []))]

    elif method == 'qpcr':
        masked = mask_sequence(tmp.sequence, tmp.mask)
        
        primers = []
        for constraints in tmp.apply('qpcr', db, params):
            # print(constraints)
            x = [p for p in next(design_primers(masked, constraints, params, []))]
            primers.extend(x)

    elif method == 'mrna':
        tmp.mrna = reconstruct_mrna(tmp.feat, db, genome, vardbs)
        constraints = tmp.apply('mrna', db, params)
        masked = tmp.mrna[0]
        primers = [p for p in next(design_primers(masked, constraints, params, []))]
    
    else:
        raise ValueError('Method is not implemented, exit.')
    
    results, aln = check_for_multiple_amplicons(primers, fp_genome)
    msg = f'{len(results)}/{len(primers)} primer pairs pass all filters'
    print(log(msg))
    st.write(msg)
    return results, tmp, aln


# https://docs.streamlit.io/knowledge-base/using-streamlit/caching-issues
# https://discuss.streamlit.io/t/unhashabletype-cannot-hash-object-of-type-thread-local/1917
@st.cache(allow_output_mutation=True)
def housekeeping(fp_config):
    print(log('Housekeeping ...'))

    with open(fp_config, 'r') as file:
        params = json.load(file)

    fp_genome = params['data']['reference']
    fp_coords = params['data']['coordinates']
    vardbs = params['data']['variation']

    # We load the annotation data later (see comment in main fn) but check
    # existance here.
    fp_annotation = params['data']['annotation']

    x = [i for i in vardbs.values()]
    for fp in [fp_coords, fp_annotation, fp_genome] + x:
        p = Path(fp)
        assert p.exists(), f'Path {fp} does not exist'

    print(log('Loading reference genome sequence'))
    genome = Fasta(fp_genome)
    print(log('Loading transcript coordinate mappings'))
    hdp = JSONDataProvider([fp_coords])

    return genome, hdp, params, vardbs


# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--config', required=True, type=str)
args = parser.parse_args()


def main(fp_config):
    st.set_page_config(layout='centered')

    genome, hdp, params, vardbs = housekeeping(fp_config)

    # Why is load annotation data here, not in housekeeping fn?
    # Cannot be opened by housekeeping bc/ second iteration will cause:
    # sqlite3.ProgrammingError: SQLite objects created in a thread can only be used in that same thread. The object was created in thread id 123145481936896 and this is thread id 123145532841984.
    db = gffutils.FeatureDB(params['data']['annotation'], keep_order=True)

    # TODO:
    # https://github.com/phiweger/primer4/issues/2
    # Why does it trigger reload in housekeeping fn but not here?
    # Something w/ session state?
    # https://stackoverflow.com/questions/64698788/streamlit-how-to-store-value-of-variable-in-cache
    # https://github.com/biocommons/biocommons.seqrepo
    x = params['data']['sequences']
    if Path(x).exists():
        print(log('Will use local transcript sequence data'))
        # TODO: It's this line:
        os.environ['HGVS_SEQREPO_DIR'] = x
    else:
        # If env var is not set, hgvs library will default to API usage, ie
        # internet connection is needed.
        print(log('Will use API to obtain sequence data'))


    st.markdown(
        r'''
        ## ❤️ Primer4U ❤️

        Example queries:

        ```bash
        # Sanger; HGVS syntax
        NM_000546.6:c.215C>G
        NM_005585.4:c.1419dup
        # qPCR; eg "::3" means we target exon 3
        NM_000546.6::4
        # mRNA; anchor primers in two exons
        NM_000546.6::5::6
        ```
        '''
        )

    # The menu
    # https://docs.streamlit.io/library/api-reference/widgets
    # Form
    # https://docs.streamlit.io/library/api-reference/control-flow/st.form
    with st.form('my_form'):
        
        order = st.text_input('Query', '')
        # Based on the selected method (PCR, ...) we'd like to update the
        # default values for the fields below. In theory, this is possible:
        # https://discuss.streamlit.io/t/circular-connection-of-slider-and-text-input/11015/4
        # BUT. Not inside a form:
        # streamlit.errors.StreamlitAPIException: With forms, callbacks can 
        # only be defined on the `st.form_submit_button`. Defining callbacks on 
        # other widgets inside a form is not allowed.
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            method = st.selectbox('Method', ('Sanger', 'qPCR', 'mRNA'))
            method = method.lower()
        with col2:
            amplicon_len_min = st.number_input('min length [bp]', value=250)
        with col3:
            amplicon_len_max = st.number_input('max length [bp]', value=600)
        with col4:
            max_variation = st.number_input('Allele frequency [%]', min_value=0., max_value=100., value=0., step=0.01) / 100
    
        # Every form must have a submit button.
        submitted = st.form_submit_button(
            'Run',
            on_click=warn(method, params, amplicon_len_min, amplicon_len_max))
        
        if submitted:
            if not order:
                st.write('Please provide a query')
                return None
            else:
                code = order.split('::')  # case: NM_000546.6::3
                # This would otherwise fail later on HGVS parsing error
                if len(code) > 1 and method == 'sanger':
                    raise ValueError('Wrong query syntax, should be something like "NM_000546.6:c.215C>G"')

                # Replace transcript version if necessary
                tx = code[0].split(':')[0]  # case : NM_005585.4:c.1419dup
                
                print(log('Sync transcript'))
                used_tx = sync_tx_with_feature_db(tx, db)
                if used_tx != tx:
                    st.warning(f'Used trancript {used_tx}')
                    code = [code[0].replace(tx, used_tx)] + code[1:]

                print(log('Primer design'))
                with st.spinner(text='In progress ...'):
                    primers, tmp, aln = gimme_some_primers(
                        method,
                        code,
                        params['data']['reference'],
                        genome,
                        hdp,
                        db,
                        vardbs,
                        params,
                        max_variation)
                st.write('Done.')

        else:
            return None
    # What's in "primers"?
    #import pdb
    #pdb.set_trace()
    # if primers: dir(primers[0])
    # dir(tmp)
    # [... 'data', 'fwd', 'insert', 'name', 'penalty', 'rev', 'save', 'to_c', 
    # 'to_g']
    # primers[0].to_g(tmp)
    # [7579197, 7579217, 7579679, 7579699]
    # TODO: coding coords

    # TODO def prepare_df() in vis
    
    df = primers_to_df(primers, tmp, order, aln)
    if df.empty:
        st.write('No primers found under the provided constrains. Relax (the constraints)!')

    else:
        # Plot something
        # data = prepare_mock_data_for_vis()
        #_ = Beauty(data).plot()
        
        # "image" is generated like so:
        # from PIL import Image
        # image = Image.open(fp)
        image = prepare_data_for_vis(tmp.data, tmp, primers)

        # Center image
        col1, col2, col3 = st.columns([1, 1000, 1])
        st.text('\n')
        with col1:
            st.write('')
        with col2:
            st.image(image)
        with col3:
            st.write('')    
        #from PIL import Image
        #image = Image.open(img_fp)
        # st.image(image)
        #img_fp.cleanup()


        # https://docs.streamlit.io/library/api-reference/data/st.dataframe
        # st.table(df)
        st.text('\n')
        st.dataframe(
            df.style.format(
                {
                    'penalty': '{:.2f}',
                    'fwd GC': '{:.2f}',
                    'rev GC': '{:.2f}',
                    'fwd Tm': '{:.2f}',
                    'rev Tm': '{:.2f}',
                }))
        # TODO: Style columns?
        # https://stackoverflow.com/questions/41654949/pandas-style-function-to-highlight-specific-columns
        # st.dataframe(df.style.highlight_max(axis=1))
        
        # Download
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-pandas-dataframe-csv
        # https://docs.streamlit.io/knowledge-base/using-streamlit/how-download-file-streamlit
        @st.cache
        def convert_df(df):
            return df.to_csv(index=False, quoting=None).encode('utf-8')

        #import pdb
        #pdb.set_trace()
        csv = convert_df(df)

        # This reloads the entire page, see "hacks.py"
        # st.download_button(
        #     "Download",
        #     csv,
        #     "file.csv",
        #     "text/csv",
        #     key='download-csv'
        #     )

        #import pandas as pd
        
        download_button_str = download_button(csv, 'primers.csv', f'Download')
        st.markdown(download_button_str, unsafe_allow_html=True)

    return None


if __name__ == '__main__':
    main(args.config)

