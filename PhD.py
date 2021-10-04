import pandas as pd
from IPython.display import display

# import all designs and standardise before use to merge with oocyst and weight tables
# selecting for columns: EH_ID, primary_infection, challenge_infection,

design_E1 = pd.read_csv(
    "https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data/Experimental_design/E57_xxxxx_Eim_DESIGN.csv"
)
design_E1.columns
