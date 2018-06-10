

def process_qPCR(qPCRfile,control_primer='LEU3',control_plasmid='CEN'):

    import pandas as pd
    import numpy as np

    #load data from file
    qPCR = pd.read_csv(qPCRfile)

    ###############average technical replicates
    Plasmid,Primer,Replicate,Ct = [],[],[],[]
    #process each clone
    for plasmid in qPCR.Plasmid.unique():
        this_plasmid = qPCR.loc[qPCR.Plasmid == plasmid]
        #process each primer for this clone
        for primer in this_plasmid.Primer.unique():
            this_primer = this_plasmid.loc[this_plasmid.Primer == primer]
            #process each biological replicate for this clone-primer combination
            #create a single line with averaged Ct values
            for replicate in this_primer.Replicate.unique():
                Plasmid.append(plasmid)
                Primer.append(primer)
                Replicate.append(replicate)
                if this_primer.loc[this_primer.Replicate == replicate].shape[0] >1:
                    Ct.append(np.mean(this_primer.loc[this_primer.Replicate == replicate].Ct))
                else:
                    Ct.append(this_primer.loc[this_primer.Replicate == replicate].Ct)

    tech_averaged = pd.DataFrame({'Plasmid':Plasmid,'Primer':Primer,'Replicate':Replicate,'Ct':Ct})

    ###############for each non-LEU3 Ct value, calculate the delta Ct with the matching LEU3 Ct value
    Plasmid,Replicate, Primer, Ct_LEU3 = [],[],[],[]
    #go through each clone
    for plasmid in tech_averaged.Plasmid.unique():
        this_plasmid = tech_averaged.loc[tech_averaged.Plasmid == plasmid]
        #go through each primer for this clone (required because multiple primers for wt)
        for primer in tech_averaged.Primer.unique():
            if primer != 'LEU3':
                this_specific_set = this_plasmid.loc[this_plasmid.Primer == primer]
                this_control_set = this_plasmid.loc[this_plasmid.Primer == control_primer]
            #go through each replicate
            for replicate in this_specific_set.Replicate.unique():
                if primer != 'LEU3':
                    Plasmid.append(plasmid)
                    Primer.append(primer)
                    Replicate.append(replicate)
                    dCt = this_control_set.loc[this_control_set.Replicate == replicate].Ct.values[0] - this_specific_set.loc[this_specific_set.Replicate == replicate].Ct.values[0]
                    Ct_LEU3.append(dCt)

    delta_Ct = pd.DataFrame({'Plasmid':Plasmid,'Primer':Primer,'Replicate':Replicate,'deltaCt_LEU3':Ct_LEU3})

    ###############calculate the ddCt value for optimised genes and the average of the corresponding wt control
    Plasmid, Replicate, deltadeltaCt_CEN = [],[],[]
    #split delta_Ct into two tables for clones and wt control
    delta_plasmids = delta_Ct.loc[delta_Ct.Plasmid != control_plasmid]
    delta_CEN = delta_Ct.loc[delta_Ct.Plasmid == control_plasmid]
    #go through each plasmid in the plasmids table
    for plasmid in delta_plasmids.Plasmid.unique():
        this_plasmid = delta_plasmids.loc[delta_plasmids.Plasmid == plasmid]
        this_CEN = delta_CEN.loc[delta_CEN.Primer == this_plasmid.Primer.iloc[0]]
        #go through each replicate for this clone
        for replicate in this_plasmid.Replicate.unique():
            Plasmid.append(plasmid)
            Replicate.append(replicate)
            ddCt = np.mean(this_CEN.deltaCt_LEU3) - this_plasmid.loc[this_plasmid.Replicate == replicate].deltaCt_LEU3.values[0]
            deltadeltaCt_CEN.append(ddCt)

    #convert ddCt to fold change opt vs wt
    fold_change_CEN = 2**-np.array(deltadeltaCt_CEN)
    return pd.DataFrame({'Plasmid':Plasmid,'Replicate':Replicate,'ddCt_CEN':deltadeltaCt_CEN, 'fold_change_CEN': fold_change_CEN})
