SELECT 
H.Survey AS Survey, 
H.HKey,
H.Year, 
H.Month, 
H.Day, 
Switch([H.Year]>2017,H.LatDegEnd+H.LatMinEnd/60,[H.Year]<2018,H.LatDegStart+H.LatMinStart/60) AS LatDeep, 
Switch([H.Year]>2017,-(H.LonDegEnd+H.LonMinEnd/60),[H.Year]<2018,-(H.LonDegStart+H.LonMinStart/60)) AS LonDeep, 
Switch([H.Year]>2017,H.LatDegStart+H.LatMinStart/60,[H.Year]<2018,H.LatDegEnd+H.LatMinEnd/60) AS LatShallow, 
Switch([H.Year]>2017,-(H.LonDegStart+H.LonMinStart/60),[H.Year]<2018,-(H.LonDegEnd+H.LonMinEnd/60)) AS LonShallow, 
H.Transect, 
H.Tr_Length AS Transect_length, 
H.QuadratSkipping, 
Q.Quadrat,
Q.Substrate1,
Q.Substrate2,
Q.CorDepthM
FROM ([1-Headers] H
	INNER JOIN [2-Quadrat] Q ON H.HKey = Q.HKey) 
ORDER BY H.Year, H.Transect, Q.Quadrat;
