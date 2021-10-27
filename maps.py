from unittest import case
import dash_bio as dashbio
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
import scipy
from scipy.stats.kde import gaussian_kde
import pandas as pd
import numpy as np
import dash_bootstrap_components as dbc


def preprocess_singlecelldata(filename, list):
    data = pd.read_csv(filename, index_col=0)
    data = data.T
    data.columns = [x.upper() for x in data.columns]
    return data[data.columns.intersection(list)]


df_alz_f = pd.read_csv('GSE44768_CR_alz_female_reduced.csv', index_col=0)
df_alz_m = pd.read_csv('GSE44768_CR_alz_male_reduced.csv', index_col=0)
df_nd_f = pd.read_csv('GSE44768_CR_nd_female_reduced.csv', index_col=0)
df_nd_m = pd.read_csv('GSE44768_CR_nd_male_reduced.csv', index_col=0)
data_colname = df_alz_f.columns
singlecell_GSM2629341_wholebrain = preprocess_singlecelldata('output_stage_2_GSM2629341_AB1442.csv', data_colname)

dict_data = {0: df_alz_f, 1: df_alz_m, 2: df_nd_f, 3: df_nd_m, 4: singlecell_GSM2629341_wholebrain};


dataset_colname = ['df_alz_f', 'df_alz_m', 'df_nd_f', 'df_nd_m', 'singlecell_GSM2629341_wholebrain']

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Histogram", href="/histogram")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("More plots...", header=True),
                dbc.DropdownMenuItem("Heatmap and clustering...", href="/heatmap"),
                # dbc.DropdownMenuItem("Histogram", href="/histogram"),
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
    ],
    brand="Data Visualization Generator",
    brand_href="#",
    color="primary",
    dark=True,
)

heatmap = html.Div([
    html.Div(dcc.Dropdown(
        id='clustergram-input-1',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(dataset_colname, dataset_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=dataset_colname[4],
        multi=False,
        placeholder='Please select dataset: ',
        style={'width': '80%'},
    ), style=dict(width='50%', display='inline-block')),
    html.Div(dcc.Dropdown(
        id='clustergram-input-2',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(dataset_colname, dataset_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=dataset_colname[4],
        multi=False,
        placeholder='Please input Gene name: ',
        style={'width': '80%'},

    ), style=dict(width='50%', display='inline-block')),

    html.Div([
        # html.P("Bandwidths: ", style={'display': 'inline-block'}),
        html.Div(id='clustergram-1'),
    ], style={'width': '50%', 'display': 'inline-block'},
    ),

    html.Div([
        # html.P("Bandwidths: ", style={'display': 'inline-block'}),
        html.Div(id='clustergram-2'),
    ], style={'width': '50%', 'display': 'inline-block'},
    ),
], style={'width': '95%'})

histogram = html.Div([
    html.Div(dcc.Dropdown(
        id='dataset',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(dataset_colname, dataset_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=dataset_colname[0],
        multi=False,
        placeholder='Please select dataset: ',
        style={'width': '80%'},
    ), style=dict(width='50%', display='inline-block')),
    html.Div(dcc.Dropdown(
        id='gene_name',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(data_colname, data_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=data_colname[1],
        multi=False,
        placeholder='Please input Gene name: ',
        style={'width': '80%'},

    ), style=dict(width='50%', display='inline-block')),
    html.Div([
        # html.P("Bin sizes: ", style={'display': 'inline-block'}),
        dcc.Graph(id="histo", clear_on_unhover=True),
        dcc.Slider(
            id='my-slider-histogram',
            min=0,
            max=0.5,
            step=0.01,
            value=0.15,
        ),
        html.Div(id='slider-output-container-histogram')
    ], style={'width': '50%', 'display': 'inline-block'},
    ),

    html.Div([
        # html.P("Bandwidths: ", style={'display': 'inline-block'}),
        dcc.Graph(id="KDE"),
        dcc.Slider(
            id='my-slider-KDE',
            min=0,
            max=0.5,
            step=0.01,
            value=0.15,
        ),
        html.Div(id='slider-output-container-KDE')
    ], style={'width': '50%', 'display': 'inline-block'},
    ),

    # second part
    html.P("Comapring..."),
    html.Div(dcc.Dropdown(
        id='dataset_2',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(dataset_colname, dataset_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=dataset_colname[0],
        multi=False,
        placeholder='Please select dataset: ',
        style={'width': '80%'},
    ), style=dict(width='50%', display='inline-block')),
    html.Div(dcc.Dropdown(
        id='gene_name_2',
        options=[{'label': s[0], 'value': str(s[1])}
                 for s in zip(data_colname, data_colname)],
        # options=[{'label': 0, 'value': str(0)}],
        value=data_colname[1],
        multi=False,
        placeholder='Please input Gene name: ',
        style={'width': '80%'},

    ), style=dict(width='50%', display='inline-block')),
    html.Div([
        # html.P("Bin sizes: ", style={'display': 'inline-block'}),
        dcc.Graph(id="histo_2", clear_on_unhover=True),
        dcc.Slider(
            id='my-slider-histogram_2',
            min=0,
            max=0.5,
            step=0.01,
            value=0.15,
        ),
        html.Div(id='slider-output-container-histogram_2')
    ], style={'width': '50%', 'display': 'inline-block'},
    ),

    html.Div([
        # html.P("Bandwidths: ", style={'display': 'inline-block'}),
        dcc.Graph(id="KDE_2"),
        dcc.Slider(
            id='my-slider-KDE_2',
            min=0,
            max=0.5,
            step=0.01,
            value=0.15,
        ),
        html.Div(id='slider-output-container-KDE_2')
    ], style={'width': '50%', 'display': 'inline-block'},
    ),

], style={'width': '95%'})

content = html.Div(id="page-content", children=[])

app.layout = html.Div([
    navbar,
    html.P("To use: simply drag the slider bar below thet graph..."),
    content,
    dcc.Location(id="url-index", refresh=False),
], style={'width': '100%'})


@app.callback(
    dash.dependencies.Output('slider-output-container-histogram', 'children'),
    [dash.dependencies.Input('my-slider-histogram', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)


@app.callback(
    dash.dependencies.Output('slider-output-container-KDE', 'children'),
    [dash.dependencies.Input('my-slider-KDE', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)


@app.callback(
    dash.dependencies.Output('slider-output-container-histogram_2', 'children'),
    [dash.dependencies.Input('my-slider-histogram_2', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)


@app.callback(
    dash.dependencies.Output('slider-output-container-KDE_2', 'children'),
    [dash.dependencies.Input('my-slider-KDE_2', 'value')])
def update_output(value):
    return 'You have selected "{}"'.format(value)


@app.callback(
    dash.dependencies.Output('histo', 'figure'),
    [dash.dependencies.Input('my-slider-histogram', 'value'), Input('dataset', 'value'), Input('gene_name', 'value')])
def update_graph_histo(bin_sizer, dataset, gene_name):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        fig = go.Figure(data=[go.Histogram(x=df[gene_name],
                                           xbins=dict(
                                               # start='1969-11-15',
                                               # end='1972-03-31',
                                               size=bin_sizer),  # 4 months bin size
                                           autobinx=False,
                                           marker_color='rgb(50, 150, 250)',
                                           )])
        # Create distplot with custom bin_size
        # fig = ff.create_distplot(hist_data, group_labels, bin_size=1, curve_type='normal',)
        return fig


@app.callback(
    dash.dependencies.Output('KDE', 'figure'),
    [dash.dependencies.Input('my-slider-KDE', 'value'), Input('dataset', 'value'), Input('gene_name', 'value')])
def update_graph_KDE(bw, dataset, gene_name):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        hist_data = [df[gene_name]]
        group_labels = [gene_name]

        fig = ff.create_distplot(hist_data, group_labels, bin_size=bw, curve_type="kde",  colors= ['rgb(50, 150, 250)']
                                 # curve_type='normal',
                                 )
        return fig


@app.callback(
    dash.dependencies.Output('histo_2', 'figure'),
    [dash.dependencies.Input('my-slider-histogram_2', 'value'), Input('dataset_2', 'value'),
     Input('gene_name_2', 'value')])
def update_graph_histo(bin_sizer, dataset, gene_name):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        fig = go.Figure(data=[go.Histogram(x=df[gene_name],
                                           xbins=dict(
                                               # start='1969-11-15',
                                               # end='1972-03-31',
                                               size=bin_sizer),  # 4 months bin size
                                           autobinx=False,
                                           marker_color='rgb(191, 52, 52)',
                                           )])
        # Create distplot with custom bin_size
        # fig = ff.create_distplot(hist_data, group_labels, bin_size=1, curve_type='normal',)
        return fig


@app.callback(
    dash.dependencies.Output('KDE_2', 'figure'),
    [dash.dependencies.Input('my-slider-KDE_2', 'value'), Input('dataset_2', 'value'), Input('gene_name_2', 'value')])
def update_graph_KDE(bw, dataset, gene_name):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        hist_data = [df[gene_name]]
        group_labels = [gene_name]

        fig = ff.create_distplot(hist_data, group_labels, bin_size=bw, curve_type="kde", colors= ['rgb(191, 52, 52)']
                                 # curve_type='normal',
                                 )
        return fig


@app.callback(
    Output('clustergram-1', 'children'),
    [Input('clustergram-input-1', 'value')])
def update_clustergram(dataset):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        # df = pd.read_csv(
        #     'https://raw.githubusercontent.com/plotly/datasets/master/Dash_Bio/Chromosomal/clustergram_brain_cancer.csv').set_index(
        #     'ID_REF')

        rows = list(df.index)
        columns = list(df.columns.values)
        if len(rows) < 2:
            return "Please select at least two rows to display."

        return dcc.Graph(figure=dashbio.Clustergram(
            data=df.loc[rows].values,
            column_labels=columns,
            row_labels=rows,
            color_threshold={
                'row': 250,
                'col': 700
            },
            hidden_labels='row',
            height=800,
            width=900
        ))


@app.callback(
    Output('clustergram-2', 'children'),
    [Input('clustergram-input-2', 'value')])
def update_clustergram(dataset):
    if not dataset:
        return {}
    else:
        df = dict_data[dataset_colname.index(dataset)]
        # df = pd.read_csv(
        #     'https://raw.githubusercontent.com/plotly/datasets/master/Dash_Bio/Chromosomal/clustergram_brain_cancer.csv').set_index(
        #     'ID_REF')

        rows = list(df.index)
        columns = list(df.columns.values)
        # print(rows)
        # print(columns)
        if len(rows) < 2:
            return "Please select at least two rows to display."

        return dcc.Graph(figure=dashbio.Clustergram(
            data=df.loc[rows].values,
            column_labels=columns,
            row_labels=rows,
            color_threshold={
                'row': 250,
                'col': 700
            },
            hidden_labels='row',
            height=800,
            width=900
        ))


@app.callback(Output("page-content", "children"), [Input("url-index", "pathname")])
def render_page_content(pathname):
    if pathname == "/":
        # return html.P("This is the content of the home page!")
        return heatmap,
    # note that the comma can let teh page ??
    elif pathname == "/histogram":
        return histogram,
    elif pathname == "/heatmap":
        return heatmap,
    # If the user tries to reach a different page, return a 404 message
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )


app.run_server(debug=True, host='0.0.0.0', port=80,
               threaded=True)
