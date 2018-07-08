# import data from Healthgrades.com to download a list of all medical oncologists
# in NYC. Save as a .pkl file.

#relevant imports
import requests, json, time
from pandas.io.json import json_normalize
import pandas as pd

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

df.to_pickle('Healthgradesdoctorlist.pkl')

#
