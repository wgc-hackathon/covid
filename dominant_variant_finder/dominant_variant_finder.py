# -*- coding: utf-8 -*-
"""
Dominant Variant Finder v0.5 (2021-04-11)
@author: maddyboo
"""

import json
import pandas as pd
import matplotlib.pyplot as plt

# %% data import
# import raw data - variant separated by country

# source data from Emma B. Hodcroft. 2021.
# "CoVariants: SARS-CoV-2 Mutations and Variants of Interest."
# https://covariants.org/
# snapshot taken 08/04/2020

# load snapshot json file
with open('210408_EUClusters_data.json') as json_file:
    raw_data = json.load(json_file)

# extract list of countries to allow looping in future
country_list = raw_data['countries'].keys()

# select a country and convert json to dataframe
country = 'United Kingdom'
country_json = json.dumps(raw_data['countries'][country])
country_data = pd.read_json(country_json, convert_dates=[
                            'week']).set_index('week')

# get list of tracked variants (columns)
variants_list = country_data.columns[1:]

# output graph of sequences vs time
fig_seq, axs_seq = plt.subplots(figsize=(20, 10))
country_data.plot.line(title=country, ylabel='sequences', ax=axs_seq).legend(
    loc='center left', bbox_to_anchor=(1.0, 0.5))
#fig_seq.savefig('sequences_'+country+'.png')

# %% variant percentages
# calculate percentage composition of each variant for that analysis 'week'
# (N.B. fortnightly reported?)

# create blank to form table from
pc_data = {}

# loop to iterate through all the variant columns, calculating percentage
# contribution of variant
for variant in variants_list:
    percentage = (country_data[variant]/country_data['total_sequences'])*100
    pc_data[variant] = percentage

# combine each column of results into a new table of variant percentage
# contribution per analysis week
pc_table = pd.concat(pc_data, axis=1)

# ouput graph of percentage contribution vs time
fig_pc, axs_pc = plt.subplots(figsize=(20, 10))
pc_table.plot.line(title=country, ylabel='percentage contribution', ax=axs_pc).legend(
    loc='center left', bbox_to_anchor=(1.0, 0.5))
#fig_pc.savefig('contrib_'+country+'.png')

# %% cumalitve totals
# calculate running total number of each variant sequence
# enables removal of false positives due to small sample dataset
cum_var_total = country_data.cumsum()

# %% identify trends
# calculate change in variant from previous analysis week
pc_change_week = pc_table.diff(periods=1)

# output graph of percentage change that week
fig_change, axs_change = plt.subplots(figsize=(20, 10))
pc_change_week.plot.line(title=country, ylabel='percentage change', ax=axs_change).legend(
    loc='center left', bbox_to_anchor=(1.0, 0.5))
#fig_change.savefig('change_'+country+'.png')

# %% variants of interest
# definitions:
# - rapid increase in proportion of cases from variant
# - above specified threshold percentage contribution
# - other rules/definitions could be added? e.g. region to region spread?

# %% rapid increase in variant
rapid_increase_triggered = {}

# loop to iterate through variants and find growth rate >5% per week
for variant in variants_list:
    growth_query = str("`"+variant+"` > 5")
    rapid_weeks = pc_change_week.query(growth_query).index
    if rapid_weeks.size > 0:
        # remove trigger dates where less than 70 sequences reported to date
        # take the first date that both criteria met
        total_query = str("`"+variant+"` > 70")
        sufficient_data = cum_var_total.query(total_query).index
        rapid_date = rapid_weeks.join(
            sufficient_data, how='inner').to_series().first('1D').index
        # don't report if both criteria not met
        if rapid_date.size > 0:
            rapid_increase_triggered[variant] = rapid_date

# combine the dates when each variant first passed rapid growth trigger
rapid_increase_table = pd.DataFrame.from_dict(rapid_increase_triggered,
                                              orient="index",
                                              columns=["Rapid Trigger At"])

# %% above specified threshold percentage
threshold_triggered = {}

# loop to iterate through variants and find first date above 20% threshold
for variant in variants_list:
    query = str("`"+variant+"` > 20")
    threshold_date = pc_table.query(query).first('1D').index
    if threshold_date.size > 0:
        threshold_triggered[variant] = threshold_date

# combine the dates when each variant first passed threshold
threshold_table = pd.DataFrame.from_dict(threshold_triggered, orient="index",
                                         columns=["Threshold Met At"])

# %% reporting
# combine and tabulate warning markers
variants_interest = rapid_increase_table.join(threshold_table,
                                              how='outer')
variants_interest["Alert Count"] = variants_interest.count(axis=1)

# output variants of interest report
print(country)
print(variants_interest.sort_values(by=['Alert Count'], ascending=False))
#variants_interest.sort_values(by=['Alert Count'], ascending=False).to_csv(
#    'interest_'+country+'.csv')
