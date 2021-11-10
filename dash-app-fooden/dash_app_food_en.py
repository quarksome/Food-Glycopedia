# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 23:37:52 2020
Python script to:
    1. access Google Sheets
    2. Convert data into pandas DF
    
@author: Nikko Bacalzo
"""

import httplib2
import os
import pandas as pd
import numpy as np
from apiclient import discovery
from google.oauth2 import service_account
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

colors = {
    'background': '#FAFAC1',
    'text': '#2160FA'
}

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
    
    
# The ID and range of a google spreadsheet.
SPREADSHEET_ID = '1thjSXAup8LVer_9JfnozMJk5JxEqEr-lhU-kaGS450s'
RANGE_NAME = 'ALL_DATA!A1:T'
credsFilename = 'client_secret.json'


def get_google_sheet():
    """
    Returns all values from the target google sheet.
    Prints values from a sample spreadsheet.
    """
    
    #setup credentials
    secret_file = os.path.join(os.getcwd(), credsFilename)
    creds = service_account.Credentials.from_service_account_file(secret_file)

    service = discovery.build('sheets', 'v4', credentials=creds)

    # Call the Sheets API
    sheet = service.spreadsheets()
    result = sheet.values().get(spreadsheetId=SPREADSHEET_ID,
                                range=RANGE_NAME).execute()
    values = result.get('values', [])

    if not values:
        print('No data found.')
    else:
        return values
    

def gsheet_to_df(values):
    """ 
    Converts Google sheet API data to Pandas DataFrame
    
    Input: Google API service.spreadsheets().get() values
    Output: Pandas DataFrame with all data from Google Sheet
    """
    header = values[0]
    rows = values[1:]
    if not rows:
        print('No data found.')
    else:
        df = pd.DataFrame(columns=header, data=rows)
    return df


def json_abs_data():
    """
    Data call to google sheets api.  
    Function will handle initial cleaning steps to reduce duplication
    during other callbacks
    
    Create a copy of the data but with relative values
    """
    gsheet = get_google_sheet()
    df = gsheet_to_df(gsheet)
    
    df.columns = ['UID',
                  'all_descriptors',
                  'food_group',
                  'water_content',
                  'Glucose',
                  'Galactose',
                  'Fructose',
                  'Xylose',
                  'Arabinose',
                  'Fucose',
                  'Rhamnose',
                  'GlcA',
                  'GalA',
                  'GlcNAc',
                  'GalNAc',
                  'Mannose',
                  'Allose',
                  'Ribose',
                  'sample_mass',
                  'anthrone'
                  ]
    
    df = df.replace('', np.nan)
    df = df.dropna(subset=['all_descriptors'])
    df['all_descriptors'] = df['all_descriptors'].str.lower()
    
    json_abs = df.to_json(orient='split')
    return json_abs

    
def json_rel_data():
    """
    Create a copy of the data but with relative values
    """
    df2 = pd.read_json(json_abs_data(), orient='split')
    
    allmonoList = ['Allose', 'Mannose', 'Fructose', 'Glucose', 'Galactose',
            'Xylose', 'Arabinose', 'Ribose',
            'Fucose','Rhamnose',
            'GlcA', 'GalA', 'GlcNAc', 'GalNAc']
    
    df2['mono_sum'] = df2[allmonoList].sum(axis=1)
    df2[allmonoList] = df2[allmonoList].div(df2['mono_sum'], axis=0)
    json_rel = df2.to_json(orient='split')
    return json_rel


def filter_plot(df, monoList, filterWord):
    """
    Function to filter dataframe according to selected mono and filter keyword
    """
    #keyword filter
    filterWord = filterWord.lower()
    temp = df[df['all_descriptors'].str.contains(filterWord)]
    
    #mono filter
    columnList = ['all_descriptors', 'UID'] + monoList
    temp = temp[columnList]
    temp_json = temp.to_json(orient='split')
    
    #plotting
    xlist = [str(x) for x in temp['UID'].tolist()]
    labels = temp['all_descriptors'].tolist()
    
    fig = go.Figure(data=[go.Bar(name = mono,
            x = temp[mono],
            y = xlist,
            orientation='h',
            hovertext = labels
            )
     for mono in monoList]
                    ) 
    
    fig.update_layout(
        height=975,
        barmode='stack',
        yaxis_type='category',
        xaxis=dict(
            title='concentration (mg/mg)',
            titlefont_size=16,
            tickfont_size=14,
            mirror = "allticks",
            side = 'top'
        ),
    )
    return temp_json, fig

   
    
allmonoList = ['Allose', 'Mannose', 'Fructose', 'Glucose', 'Galactose',
            'Xylose', 'Arabinose', 'Ribose',
            'Fucose','Rhamnose',
            'GlcA', 'GalA', 'GlcNAc', 'GalNAc']


#LAYOUT OF THE DASHBOARD
app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(children='Lebrilla Lab - Food Glycomics Encyclopedia',
            style={
                'textAlign': 'center',
                'color': colors['text']
        }
    ),
    
    html.Br(),
    
    html.Label('Choose monosaccharide quantitation data:',
    ),
    
    dcc.RadioItems(
        id='quant-option',
        options=[
        {'label': 'absolute', 'value': 'absolute'},
        {'label': 'relative', 'value': 'relative'}
        ],
        value='absolute',
        labelStyle={'display': 'inline-block'}
    ),
    
    html.Br(),
    
    html.Label('Choose monosaccharides to display:',
    ),
    
    dcc.Checklist(
        id='mono-filter',
        options=[{'label': x, 'value': x} for x in allmonoList],
        value=['Glucose'],
        labelStyle={'display': 'inline-block'}
    ),
    
    html.Br(),
    
    html.Label('Enter filter keyword',
    ),

    dcc.Input(
        id='filter-word',
        placeholder='Enter a filter keyword...',
        type='text',
        value='apple',
    ),
    
    html.Br(),
    
    html.Label('Hover mouse pointer to graph for sample information.',
    ),
    
    html.Br(),
    
    html.Label('Click on graph legend to toggle display.',
    ),

    dcc.Graph(id='filtered-data-GRAPH'),
    
    html.Div(id='main-data-holder',
                 children=json_abs_data(),
                 style={'display': 'none'}
    ),
    
    html.Div(id='relative-data-holder',
             children=json_rel_data(),
             style={'display': 'none'}
    ),

    html.Div(id='filtered-data-holder', style={'display': 'none'}
    ),
    
    html.Div(id='filtered-rel-data-holder', style={'display': 'none'}
    )

    
])


@app.callback(
    [Output('filtered-data-holder', 'children'),
     Output('filtered-rel-data-holder', 'children'),
     Output('filtered-data-GRAPH', 'figure')],
    [Input('main-data-holder', 'children'),
     Input('relative-data-holder', 'children'),
     Input('quant-option', 'value'),
     Input('mono-filter', 'value'),
     Input('filter-word', 'value')])
def data_process(absData, relData, quant_option, monoList, filterWord):
    """
    This function will use the main hidden data component as input.
    Then filter the data from the specificed MONO_CHECKBOX and FILTER_KEYWORD
    
    Returns the new data in string format stored to another hidden
    to be reused in multiple locations. 
    Returns plotly figure object.
    """
    if quant_option=='absolute':
        df = pd.read_json(absData, orient='split')
        temp_abs_json, fig = filter_plot(df, monoList, filterWord)
        temp_rel_json = pd.DataFrame().to_json(orient='split')
        fig.update_layout(
            xaxis=dict(
                title='concentration (mg/mg)')
        )
    
    elif quant_option=='relative':
        df = pd.read_json(relData, orient='split')
        temp_rel_json, fig = filter_plot(df, monoList, filterWord)
        temp_abs_json = pd.DataFrame().to_json(orient='split')
        fig.update_layout(
            xaxis=dict(
                title='fractional relative abundance')
        )
        
    return temp_abs_json, temp_rel_json, fig





if __name__ == '__main__':
    app.run_server(debug=True)
