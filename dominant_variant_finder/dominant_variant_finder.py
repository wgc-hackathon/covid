# -*- coding: utf-8 -*-
"""
Dominant Variant Finder v0.1.1 (2021-04-06)
@author: maddyboo
"""

import pandas as pd

# %% data import
# import raw daily data
data_csv = 'Imaginary_Data.csv'
raw_data = pd.read_csv(data_csv, parse_dates=(['Date']), index_col=('Date'),
                       dayfirst=(True))

# get list of variants (columns)
variants_list = raw_data.columns

# back-fill period before variant existed with zero cases
backfill_data = raw_data.fillna(0)

# get total cases reported that day
backfill_data['Total'] = backfill_data.sum(axis=1)

# calculate centered 7-day rolling average to smooth data
cleaned_data = backfill_data.rolling(7, center=True).mean()
cleaned_data.plot.line()

# %% daily variant percentages
# calculate percentage composition of each variant for that day

# create blank to form table from
pc_data = {}

# loop to iterate through all the variant columns, calculating percentage
# contribution of variant to the daily new cases total
for variant in variants_list:
    percentage = (cleaned_data[variant]/cleaned_data['Total'])*100
    pc_data[variant] = percentage

# combine each column of results into a new table of variant percentage
# contribution per day
pc_table = pd.concat(pc_data, axis=1)
pc_table.plot.line()

# %% identify trends
# calculate change in variant percentage over last 7 days
pc_change_week = pc_table.diff(periods=7)
pc_change_week.plot.line()

# %% variants of interest
# definitions:
# - rapid increase in proportion of cases from variant
# - above specified threshold percentage of daily reported cases
# - other rules/definitions could be added? e.g. region to region spread?

# %% rapid increase in variant
rapid_increase_triggered = {}

# loop to iterate through variants and find growth rate >5% per week
for variant in variants_list:
    query = str(variant+" > 5")
    rapid_date = pc_change_week.query(query).first('1D').index
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
    query = str(variant+" > 20")
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
print(variants_interest.sort_values(by=['Alert Count'], ascending=False))
