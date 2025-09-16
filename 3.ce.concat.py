# %%
import pandas as pd
import os 
import argparse
from datetime import date
# %%
parser = argparse.ArgumentParser(description='concatinate all the library sequences splice detections from a single experiment.')
parser.add_argument('-p', '--pwrdir', type=str, help='path to put all the putouts')
parser.add_argument('-n', '--name', type=str, help='name of experiment - prefix')
args = parser.parse_args()
ttoday = date.today().strftime("%Y-%m-%d")
# star = pd.concat([pd.read_csv('csvfiles/' + file) for file in os.listdir('csvfiles/') if file.endswith('.csv')]).drop(columns='isoform')
pd.concat([pd.read_csv(args.pwrdir + '/csvfiles/' + file) for file in os.listdir(args.pwrdir + '/csvfiles/') if file.endswith('.csv')]).to_csv(args.name + '.' + ttoday + '.csv', index=False, header=True)
# pd.concat([pd.read_csv(args.pwrdir + '/rawsplicedata/' + file) for file in os.listdir(args.pwrdir + '/rawsplicedata/') if file.endswith('.csv')]).to_csv(args.name + '.star.raw.csv', index=False, header=True)

