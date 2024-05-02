from Tissue_obj import Tissue_obj
from Artifact_detect import Artifact_detect
import pandas as pd
import numpy as np
import copy

########
# This class contain functions removing function
#
#
#
########

class Artifact_remove(Artifact_detect):
    def __init__(self, dir):
        super().__init__(dir)
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        self.tissue_depth = 0
        self.tissue_depth_df = self.fn_spot_dis_to_edge()
        self.remove_procedure = []
        self.removed_spots = {}

        self.remaining_depth = self.tissue_depth
        

    def fn_spot_dis_to_edge(self):
        #print(0)
        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])
        df['spot_depth'] = [0] * (df.shape[0])
        # this function returns the maximum depth and a dataframe of spots & its depth level
        #print("1")

        # 1. Firstly generate a set of spots which are "covered" "on edge" or "on border".

        # Get edge spots
        tissue_position = self.tissue_position
        zero_idx = tissue_position[tissue_position.in_tissue == 0].index.to_list()
        covered_idx = tissue_position[tissue_position.in_tissue == 1].index.to_list()
        neighbors = set()
        for i in range(len(zero_idx)):
                barcode = self.tissue_position.at[zero_idx[i],'barcode']
                neighbors.update(self.get_neighbour_1_bar(barcode))
        covered_tissue = set(self.tissue_position.loc[covered_idx,'barcode'])
        #print(len(neighbors))
        # Get border spots
        df_border_spot = self.get_border()
        #print(len(df_border_spot.barcode))
        neighbors.update(df_border_spot.barcode)

        neighbors.intersection_update(covered_tissue)
        neighbors.discard(None)
        for barcode in neighbors:
            df.loc[df.barcode == barcode, "spot_depth"] = 1
        remained = covered_tissue.difference(neighbors)
        #print("2")


        # 2. starting from spot with distance 1 to the "edge" or "border"
        level = 2 
        while len(remained) != 0:
            #print(len(remained))
            #print(len(neighbors))
            ngh_iter = copy.deepcopy(neighbors)
            for barcode in ngh_iter:
                neighbors.update(self.get_neighbour_1_bar(barcode))
                
            neighbors.intersection_update(covered_tissue)
            neighbors.discard(None)
            added = neighbors.difference(ngh_iter)    
            for barcode in added:
                df.loc[df.barcode == barcode, "spot_depth"] = level
            level = level + 1
            remained = covered_tissue.difference(neighbors)
        self.tissue_depth = level 
        #print(3)
        print(f"The depth of tissue is {level}")

        return df

    def fn_spot_inclusion_condition(self):
        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])

        return df




    def remove_edge(self, distance = 1):
        self.remove_procedure.append("edge")
        # return warn "all spots are removed"
        removed_spot_n = sum( (self.tissue_depth_df.spot_depth <= distance) & (self.spot_inclusion_condition.in_tissue.astype('bool')))
        self.spot_inclusion_condition.loc[self.tissue_depth_df.spot_depth <= distance, 'in_tissue'] = 0
        print(f"We removed {removed_spot_n} edge spots with distance to edge {distance} and less")
        return 

    def remove_border(self):
        self.remove_procedure.append("border")
        df_border_spot = self.get_border()
        for barcode in df_border_spot.barcode:
            self.spot_inclusion_condition.loc[ self.spot_inclusion_condition.barcode == barcode, 'in_tisse'] = 0
        print(f"We removed {len(df_border_spot.barcode)} border spots")
        return

    def remove_malfunction(self):
        # Gel all malfunction points, add to remove_procedure list 
        self.remove_procedure.append("malfunction")
        df_outlier_spot = Artifact_detect.outlier(self.get_sum())
        for barcode in df_outlier_spot.barcode:
            self.spot_inclusion_condition.loc[ self.spot_inclusion_condition.barcode == barcode, 'in_tisse'] = 0
        n_outlier = len(df_outlier_spot.barcode)
        print(f"We removed {n_outlier} outlier spots")
        return
    



#----------- 
# Functions for check object current conditions (in our case check what spots removed, 
# procedure done for removing)

    def simple_cleanser(self): 
        self.remove_procedure = []
        self.removed_spots = {}
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        print("Back to the initial Tissue Sample")
        return 
    
    def review_removing(self):
        [print(i) for i in self.remove_procedure]
        return 
    
    def save(self, dir):
        # Save the position list dataframe 
        df = copy.deepcopy(self.tissue_position)

        df.in_tissue = self.spot_inclusion_condition.in_tissue
        df.to_csv(dir + '/tissue_position.csv', index = False)

        return