import streamlit as st
import pandas as pd
from pydeseq2.preprocessing import deseq2_norm
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

count_data = pd.read_csv('E-GEOD-22260-raw-counts.tsv', sep='\t')
metadata = pd.read_csv('E-GEOD-22260-experiment-design.tsv', sep='\t')

st.write("""
    Using a datasaet from this research [insert research] I want to identify the ten genes that show the most differences
    between a healthy individual and cancer patients. First, we load in and inventorize the data
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

metadata = metadata.set_index('Run')

metadata = metadata.drop(columns=['Sample Characteristic Ontology Term[disease]', 'Sample Characteristic[organism]',
                                  'Sample Characteristic Ontology Term[organism]', 'Sample Characteristic[organism part]',
                                  'Sample Characteristic Ontology Term[organism part]', 'Sample Characteristic Ontology Term[tumor stage]',
                                  'Factor Value[disease]', 'Factor Value Ontology Term[disease]', 'Analysed'])


count_data = count_data.set_index('Gene ID')
count_data = count_data.drop(columns=['Gene Name'])

count_entries_before = len(count_data)

count_data = count_data.loc[(count_data!=0).any(axis=1)]

count_entries_after = len(count_data)

st.write(f"""
    Dropping the non-expressed/detected genes brings the dataframe down from **{count_entries_before}** to 
    **{count_entries_after}**. The next step is data normalisation and scaling. This is to compensate for library depth 
    (which can cause difference between gene expression based on technical variables rather than technological ones) 
    and scaling to place our data into a smaller range and with it make differences more apparent.\n
    """)

### DeSeq2 requires the columns to be the genes and the rows the subjects.

count_data_transposed = count_data.transpose()
count_data_normalised_tuple = deseq2_norm(count_data_transposed)
count_data_normalised = count_data_normalised_tuple[0]
scaler = StandardScaler()
count_data_scaled = scaler.fit_transform(count_data_normalised)
count_data_scaled = pd.DataFrame(count_data_scaled, columns=count_data_normalised.columns, index=count_data_normalised.index)

count_data_tranposed_meta = count_data_transposed.join(metadata)
count_data_scaled_meta = count_data_scaled.join(metadata)

### The boxplots have been left out as they are slow to load with the data required
xvar = st.selectbox('Select x-axis:', count_data_transposed.columns)
yvar = st.selectbox('Select y-axis:', count_data_transposed.columns)

col1, col2 = st.columns(2)
with col1:
    # st.write(px.box(count_data_transposed))
    st.write(px.scatter(count_data_tranposed_meta,  x=xvar, y=yvar,
                        color='Sample Characteristic[disease]', hover_data= {'Subject' : count_data_tranposed_meta.index}))
with col2:
    # st.write(px.box(count_data_normalised))
    st.write(px.scatter(count_data_scaled_meta, x=xvar, y=yvar,
                        color='Sample Characteristic[disease]', hover_data= {'Subject': count_data_scaled_meta.index}))

st.write("""
    We can see that among the two selected genes we are now working along a smaller distribution. This should make the significant
    differences between groups easier to identify. As the next step, I want to take my normalised data and process it through a PCA to see
    if our two groups show distinction.
    """)

count_pca = PCA()
count_transformed = count_pca.fit_transform(count_data_scaled)

col_names = [f'component {i+1}' for i in range(count_transformed.shape[1])]

count_transformed_df = pd.DataFrame(count_transformed, columns=col_names, index=count_data_scaled.index)

### Need to re-add the metadata
### Ruhroh, metadata has an index that the other doesn't have
count_transformed_df = pd.concat([count_transformed_df, metadata['Sample Characteristic[disease]']], axis=1)

xvar = st.selectbox('Select x-axis:', count_transformed_df.columns)
yvar = st.selectbox('Select y-axis:', count_transformed_df.columns)

st.write(px.scatter(count_transformed_df, x=xvar, y=yvar,
                    color='Sample Characteristic[disease]', hover_data= {'Subject': count_transformed_df.index}))

st.write("""
    Among the variuos mixes of PCAs, one thing stands out. There is usually one subject that darts away from the rest, 
    but notably: this subject seems to differ between combinations, not one individual who is an outlier each time.\n
    While the distance between the two groups isn't big, they do usually seem to be moving two different directions, thus 
    there should be rna expression that could explain this.\n
    That means the next step is looking at the level of difference between the two groups. For this I will be looking at the
    log fold change.
    """)
