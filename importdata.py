# import data from Healthgrades.com to download a list of all medical oncologists
# in NYC. Save as a .pkl file.

#relevant imports
import requests, json, time
from Bio import Entrez
from pandas.io.json import json_normalize
import pandas as pd
import pickle

#import healthgrades NYC data
base_api = 'https://www.healthgrades.com/api3/usearch'
params = {
    'userLocalTime': '13:35',
    'practicing_specialties': 'ps456',
    'what': 'Oncology',
    'entityCode': 'PS592',
    'searchType': 'PracticingSpecialty',
    'spec': 67,
    'category': 'provider',
    'where': 'New York, NY',
    'pt': '40.71455, -74.007118',
    'sort.provider': 'bestmatch',
    # sessionId: S3c20f1ef717e279c
    # requestId: Rdbfc23779ea047aa
    'pageNum': 1,
    'isFirstRequest': 'true',
    'debug': 'false',
    'isAtlas': 'false'
}

itemsPerPage = 36
resp = requests.get(base_api, params=params)
resp = json.loads(resp.text)
numPages = int(resp['search']['searchResults']['totalCount'] / itemsPerPage) + 1
for page in range(1, numPages + 1):
    params['pageNum'] = page
    resp = requests.get(base_api, params=params)
    time.sleep(0.1)
    resp = json.loads(resp.text)
    providers = resp['search']['searchResults']['provider']
    data = (providers['results'])
    if not 'df' in locals():
        df=pd.DataFrame.from_dict(json_normalize(data),orient='columns')
    else:
        df_temp=pd.DataFrame.from_dict(json_normalize(data),orient='columns')
        df = pd.concat([df,df_temp])

#save file
df.to_pickle('Healthgradesdoctorlist.pkl')

#add imports from Pubmed and Scopus
#ccreate a dictionary containing pubmed IDs of all articles published by doctors
pub_dict = {}
#read email address
file_object = open('email.key', 'r')
email = file_object.read()
file_object.close()
Entrez.email = file_object
count = 0
for name in df.head(1)['fullname']:
    if not name in publist.keys():
        searchterm = '(cancer[MeSH Terms]) AND ' + name + '[Author] '
        handle = Entrez.esearch(db="Pubmed", retmax=500, term= searchterm, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        pub_dict[name] = record['IdList']
        time.sleep(0.4)
        count+=1
#save as a .pkl file
output = open('pubmeddict.pkl', 'wb')
pickle.dump(publist, output)

#make a big list of all pubmed IDs (PMID)
#get article list and combine into big list
uidlist_master = []
for name in pub_dict:
    uidlist_master += pub_dict[name]


#Get full abstracts for all pubmed IDs (PMID) in uidlist_master
handle = Entrez.efetch(db='pubmed', id=uidlist_master, retmode='xml', rettype='abstract')
record = Entrez.read(handle)
handle.close()


extract articletitle, journal title, year and abstract
#Use dictionaries to be very careful
articletitle = {}
journaltitle = {}
journalyear = {}
abstract = {}

for article in record['PubmedArticle']:
    ID = (article['MedlineCitation']['PMID'])
    articletitle[ID] = article['MedlineCitation']['Article']['ArticleTitle']
    journaltitle[ID] = (article['MedlineCitation']['Article']['Journal']['Title'])
    if 'Year' in article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
        journalyear[ID] = (article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'])
    if 'Abstract' in article['MedlineCitation']['Article']:
        abstract[ID] = (article['MedlineCitation']['Article']['Abstract']['AbstractText'][0])

#create a dataframe of with four columns for each publication ID
df_temp1 = pd.DataFrame.from_dict(articletitle,orient = 'index', columns = ['ArticleTitle'])
df_temp2 = pd.DataFrame.from_dict(journaltitle,orient = 'index', columns = ['JournalTitle'])
df_temp3 = pd.DataFrame.from_dict(journalyear,orient = 'index', columns = ['JournalYear'])
df_temp4 = pd.DataFrame.from_dict(abstract,orient = 'index', columns = ['Abstract'])
pub_df = pd.concat([df_temp1,df_temp2,df_temp3,df_temp4], axis=1, sort=False)



#Use SCOPUS API to retreive citations for each pubmed ID
file_object = open('scopusapi.key', 'r')
api_key = file_object.read()
file_object.close()

#create a dictionary to hold the citations
citations = {}
base_api = 'https://api.elsevier.com/content/search/scopus'
for PMID in uid_df.index:
    query_param = 'PMID(' + PMID + ')'
    params = {
        'apiKey': api_key,
        'query': query_param,
        'field': 'citedby-count'
    }
    response = requests.get(base_api, params=params)
    response = json.loads(response.text)
    if 'citedby-count' in response['search-results']['entry'][0].keys():
        citations[PMID] = float(response['search-results']['entry'][0]['citedby-count'])
    time.sleep(0.1)

df_temp1 = pd.DataFrame.from_dict(citations,orient = 'index')
df_temp2 = df_temp1.rename(index=str, columns={0: "Citation Number"})
pub_df = pd.concat([pub_df,df_temp2], axis=1)
#save as a .pkl file
pub_df.to_pickle('ArticleDetails.pkl')
