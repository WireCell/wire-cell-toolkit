import torch
from roi_metrics import roi_table, roi_efficiency, roi_purity                                                                                                                           
lab = torch.zeros(1, 2, 8)                                                                                                                                                              
lab[0,0,1:3] = 1   # true ROI A                                                                                                                                                         
lab[0,0,5:7] = 1   # true ROI B                                                                                                                                                         
lab[0,1,7]   = 1   # true ROI C                                                                                                                                                         
y = torch.zeros(1, 2, 8)                                                                                                                                                                
y[0,0,2] = 0.9     # inside A                                                                                                                                                           
y[0,0,4] = 0.9     # isolated noise ROI (gap at 3 keeps it separate)                                                                                                                    
y[0,0,7] = 0.9     # end of row 0                                                                                                                                                       
y[0,1,0] = 0.9     # start of row 1 -> separate ROI, no wrap                                                                                                                            
y[0,1,7] = 0.6     # inside C
r = roi_table(y>0.5, lab==1)
print("reco ROIs:", {k: v.tolist() for k, v in r.items()})
print("eff", float(roi_efficiency(y, lab)), "expect", round(2/3,4))
print("pur", float(roi_purity(y, lab)),     "expect", 2/5)

print('true:\n', lab==1)
print('reco:\n', y>0.5)
