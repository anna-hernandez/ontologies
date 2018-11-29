import pandas as pd
import glob
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import json
from pandas.io.json import json_normalize
import requests
import json
import time
import re

from bokeh.io import show, output_file, output_notebook
from bokeh.plotting import figure
from bokeh.models.graphs import from_networkx
from bokeh.models import tools, Circle, Plot, Range1d, HoverTool
from bokeh.palettes import Viridis11


# creates psimi ontology terms
class Ontology:
    
    """Creates instances of the Ontology class.

    For an input ontology code, it recursively finds all children nodes up to the 
    specified depth level in the ontology.

    Requires Pandas, Numpy, Json, Requests, and Time.

    Parameters
    ----------
    params : dictionary of parameters (parameter_name:value)
        The dictionary will be stored as an instance property and parameters used in methods.

    url : ************** UPDATE ********************

    ontology_id : string
        Acronym of the queried ontology.

    children_levels : int or None (default=None)
        Levels to query starting at the current code level (level 0). Level 1 corresponds to children,
        2 to grand-children, etc.

    Returns
    -------
    df : dict
        A dictionary of positions keyed by node

    Examples
    --------
    >>> 

    """
    def __init__(self, params):
        self.request_date = ''
        self.last_update = ''
        self.name = ''
        self.df = pd.DataFrame()
        self.params = params

    def __str__(self):
        return 'Instance of "Ontology" class'
    
    
    
    
    # ***********************
    # collect descriptive data
    def info(self):
        url = 'https://www.ebi.ac.uk/ols/api/ontologies/{}'.format(self.params['ontology_id'])
        try:
            r = requests.get(url)
            if r.status_code == 200:
                self.request_date = r.headers['Date']
                r_json = json.loads(r.text)
                self.last_update = r_json['updated']
                self.name = r_json['config']['title']
            else:
                print('Status:', r.status_code)
        except Exception as e:
            print(e)
    
    
    
    # ***********************
    def pull_hierarchy(self):
        
        # inner functions
        # ---------------
        def get_codes(self): 

            # collect psimi terms
            ids, labels, definitions, \
            next_children, parents, new_parents = [np.empty((0,0)),np.empty((0,0)), np.empty((0,0)),\
                                                   np.empty((0,0)),np.empty((0,0)),np.empty((0,0))]
            
            for url in self.params['url']:
                try:

                    r2 = requests.get('https://www.ebi.ac.uk/ols/api/ontologies/{}/'.format(self.params['ontology_id']) + url[1])
                    if r2.status_code == 200:

                        r2_json = json.loads(r2.text)

                        ontology_terms = r2_json['_embedded']['terms']


                        for term in ontology_terms:

                            try:
                                term_id = term['annotation']['id'][0]
                                ids = np.append(ids, term_id)
                                if term_id == url[0]:
                                    parents = np.append(parents, 'root')
                                else:
                                    parents = np.append(parents, url[0])
                            except KeyError:
                                continue

                            labels = np.append(labels, term['label'])
                            if 'definition' in term['annotation'].keys():
                                term_definition = term['annotation']['definition'][0].replace('\n',' ')
                                definitions = np.append(definitions, term_definition)
                            else:
                                definitions = np.append(definitions, np.nan)

                            # store link to next level of children terms
                            # at the root level there will be one single link to store. As you go down the 
                            # hierarchical graph, more links will appear, up to one for each new term currently explored
                            # this is a temporary array, so I don't store it as class property
                            if 'hierarchicalChildren' in term['_links'].keys():
                                href = term['_links']['hierarchicalChildren']['href']
                                href = '/'.join(href.split('/')[7:])
                                next_children = np.append(next_children, href)
                                new_parents = np.append(new_parents, term['annotation']['id'][0])

                        zipped_urls = list(zip(new_parents, next_children))        

                except requests.ConnectionError:
                    print("Connection refused. Waiting 5 seconds before trying again.")
                    time.sleep(5)
                    break
                    print('Status:', r1.status_code)

            # prepare dataframe
            zipped_urls = list(zip(new_parents, next_children))        
            codes_df = pd.DataFrame(list(zip(ids, labels, definitions, parents)), columns = ['code','label','description', 'parent'])
            if next_children is not None:
                return codes_df, zipped_urls
            else:
                return codes_df
        # ---------------
        
        # ---------------
        def find_children(self):
            
            current_level_data = get_codes(self)

            if current_level_data is not None:

                self.df = pd.concat([self.df, current_level_data[0]])
                if len(current_level_data[1]) > 0:

                    # update params property
                    self.params['url'] = current_level_data[1]
                    return True
                else:
                    return False
        # ---------------
        # end inner functions

        # outter function main body
        flag = True
        if self.params['children_levels'] == 'all':
            while flag == True:
                flag = find_children(self)
                
        else:
            level_tracker = 0
            while level_tracker < params['children_levels']:
                flag = find_children(self)
                level_tracker +=1

                
        # reset index from 1 to insert into SQL dataframe
        self.df.index = np.arange(1, self.df.shape[0]+1)
        
    def calculate_depth(self):
        
        def check_parent(aString):
            if aString == 'root':
                return True
            else:
                return False

        levels = np.empty((0,0))
        for i,r in self.df.iterrows():
            parent = r.parent
            lev = 0
            while check_parent(parent) == False:
                parent = self.df[self.df.code == parent].parent.values[0]
                lev += 1
            levels = np.append(levels, lev)

        self.df['depth'] = levels
        
        
    def plot(self): 
        
        # requires depth column
        if 'depth' not in self.df.columns:
            self.calculate_depth()
        else:
            G = nx.Graph()
            G.add_edges_from([tuple(x) for x in self.df[['code','parent']].values if x[1] != 'root'])
            G.add_nodes_from([x for x in  np.unique(self.df[['code', 'parent']].values) if x != 'root'])
            labels= {x[0]:x[1] for x in  self.df[['code', 'label']].values}
            depth = {x[0]:x[1] for x in  self.df[['code', 'depth']].values}
            codes = {x[0]:x[1] for x in  self.df[['code', 'code']].values}
            nx.set_node_attributes(G, depth, 'depth')
            nx.set_node_attributes(G, codes, 'code')
            nx.set_node_attributes(G, labels, 'label')
            #colors = {x[0]:Viridis11[int(x[1])] for x in  psimi.df[['code', 'depth']].values}
            colors = [Viridis11[int(x[1])] for x in  self.df[['code', 'depth']].values]
            nx.set_node_attributes(G, colors, 'color')
        
            # create figure panel
            fig = figure(title = 'Detection Methods Ontology',
                         plot_width = 600,
                         plot_height = 600,
                         x_range=(-3,3), y_range=(-3,3))


            options = {
                'node_color': 'red',
                'node_size': 100,
                'width': 3}

            # create graph layout
            graph_renderer = from_networkx(G, nx.kamada_kawai_layout, scale=3, center=(0,0))

            # append graph to figure
            fig.renderers.append(graph_renderer)


            output_file('test.html')
            show(fig)


# ----------------------------- MAIN ----------------------------- #

psimi = Ontology(params={'url':[('MI:0000','/terms?iri=http://purl.obolibrary.org/obo/MI_0000')],
                        'ontology_id':'mi',
                        'children_levels':'all'})
psimi.info()
psimi.pull_hierarchy()
psimi.calculate_depth()
psimi.plot()

psimi.df.to_csv(psimi.name,sep='\t')


