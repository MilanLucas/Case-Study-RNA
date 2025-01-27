import streamlit as st
import pandas as pd

count_data = pd.read_csv('E-GEOD-22260-raw-counts.tsv', sep='\t')
metadata = pd.read_csv('E-GEOD-22260-experiment-design.tsv', sep='\t')

st.write("""
    I will be analysing a cancer data set and compare it to healthy samples. First, we load in and inventorize the data
    to see what we're working with.
    """)

col1, col2 = st.columns(2)
with col1:
    st.dataframe(count_data)
with col2:
    st.dataframe(metadata)

st.write(f"""
    Both dataframes are quite packed, and contain data not necessary for what we are doing. For the metadata we are only
    interested in three things: the run, whether the patient is sick, and what stage it is in. For the genes, firstly I
    will turn the Gene ID into the index and drop the Gene Names (since not all genes are named). This change will make
    it easier to filter out entries that have no expression among samples.""")

count_data = count_data.set_index('Gene ID')
count_data = count_data.drop(columns=['Gene Name'])

count_entries_before = len(count_data)

count_data = count_data.loc[(count_data!=0).any(axis=1)]

count_entries_after = len(count_data)

st.write(f"""
    Dropping the non-expressed/detected genes brings the dataframe down from **{count_entries_before}** to 
    **{count_entries_after}**.
    """)