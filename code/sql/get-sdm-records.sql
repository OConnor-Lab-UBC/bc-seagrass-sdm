SELECT 
H.Survey AS Survey, 
H.Region,
H.HKey,
H.Year, 
H.Month, 
H.Day, 
H.LatDegStart + H.LatMinStart / 60 AS LatDeep,
-(H.LonDegStart + H.LonMinStart / 60) AS LonDeep,
H.LatDegEnd + H.LatMinEnd / 60 AS LatShallow,
-(H.LonDegEnd + H.LonMinEnd / 60) AS LonShallow,
H.Transect,
H.Visibility,
H.[High vegetation],
H.RA_RSU,
H.RA_GSU,
H.RA_PSU,
RA_RSC,
PC_ZM,
PC_PH,
PC_PT
FROM [1 - Header] H;


