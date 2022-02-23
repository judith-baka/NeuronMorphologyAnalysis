import gspread
from oauth2client.service_account import ServiceAccountCredentials
from googleapiclient.discovery import build    
import pandas as pd
#import matplotlib.pyplot as plt
#import numpy as np
import os
import numpy as np# use creds to create a client to interact with the Google Drive API



scope = ['https://www.googleapis.com/auth/analytics.readonly',
      'https://www.googleapis.com/auth/drive',
      'https://www.googleapis.com/auth/spreadsheets',
      ]#['https://spreadsheets.google.com/feeds']
#creds = ServiceAccountCredentials.from_json_keyfile_name('client_secret.json', scope)
creds = ServiceAccountCredentials.from_json_keyfile_name('./credentials/online_notebook_drive_api.json', scope)
client = gspread.authorize(creds)
creds_drive = ServiceAccountCredentials.from_json_keyfile_name('./credentials/online_notebook_drive_api.json', scope)

#%% open 

def fetch_lastmodify_time(spreadsheetname):
    modifiedtime = None
    ID = None
    service = build('drive', 'v3', credentials=creds)
    wb = client.open(spreadsheetname)
    ID = wb.id
    if ID:
        modifiedtime = service.files().get(fileId = ID,fields = 'modifiedTime').execute()
    return modifiedtime

def fetch_sheet_titles(spreadsheetname):
    wb = client.open(spreadsheetname)
    sheetnames = list()
    worksheets = wb.worksheets()
    for sheet in worksheets:
        sheetnames.append(sheet.title)
    return sheetnames

def fetch_sheet(spreadsheet_name,sheet_title):
    #%%
    wb = client.open(spreadsheet_name)
    sheetnames = list()
    worksheets = wb.worksheets()
    for sheet in worksheets:
        sheetnames.append(sheet.title)
        #%%
    if sheet_title in sheetnames:
        print(sheet_title)
        idx_now = sheetnames.index(sheet_title)
        if idx_now > -1:
            params = {'majorDimension':'ROWS'}
            temp = wb.values_get(sheet_title+'!A1:OO10000',params)
            temp = temp['values']
            header = temp.pop(0)
            data = list()
            for row in temp:
                data.append(row)
            df = pd.DataFrame(data, columns = header)
            #%%
            return df
        else:
            return None
    else:
        return None

def update_metadata(notebook_name,metadata_dir):
# =============================================================================
#     lastmodify = fetch_lastmodify_time(notebook_name)
#     with open(os.path.join(metadata_dir,'last_modify_time.json')) as timedata:
#         lastmodify_prev = json.loads(timedata.read())
#     if lastmodify != lastmodify_prev:
# =============================================================================
    print('updating metadata from google drive')
    sessions = fetch_sheet_titles(notebook_name)
    for session in sessions:
        df_wr = fetch_sheet(notebook_name,session)
        if type(df_wr) == pd.DataFrame:
            df_wr.to_csv(os.path.join(metadata_dir,'{}.csv'.format(session))) 
# =============================================================================
#             #%
#     with open(os.path.join(metadata_dir,'last_modify_time.json'), "w") as write_file:
#         json.dump(lastmodify, write_file)
#         #%
#     print('metadata updated')
#     else:
#         print('metadata is already up to date')
# =============================================================================

