SELECT 
H.Survey AS Survey, 
H.HKey,
H.Year, 
H.Month, 
H.Day, 
H.LatDegStart + H.LatMinStart / 60 AS LatDeep,
-(H.LonDegStart + H.LonMinStart / 60) AS LonDeep,
H.LatDegEnd + H.LatMinEnd / 60 AS LatShallow,
-(H.LonDegEnd + H.LonMinEnd / 60) AS LonShallow,
H.Transect
FROM [1 - Header] H;

