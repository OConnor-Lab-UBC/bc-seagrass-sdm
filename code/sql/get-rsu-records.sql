SELECT
H.[Key] AS HKey,
H.Source AS Survey,
H.Location,
H.Year AS Year,
H.Month AS Month,
H.Day AS Day,
H.StatArea AS Major_stat_area_code,
H.SubArea AS Stat_subarea_code,
CASE
    WHEN H.TransOrientation = 'D' THEN H.LatDegStart + H.LatMinStart / 60
    WHEN H.TransOrientation = 'S' THEN H.LatDegEnd + H.LatMinEnd / 60
    END AS LatDeep,
CASE
    WHEN H.TransOrientation = 'D' THEN -(H.LongDegStart + H.LongMinStart / 60)
    WHEN H.TransOrientation = 'S' THEN -(H.LongDegEnd + H.LongMinEnd / 60)
    END AS LonDeep,
CASE 
    WHEN H.TransOrientation = 'D' THEN H.LatDegEnd + H.LatMinEnd / 60
    WHEN H.TransOrientation = 'S' THEN H.LatDegStart + H.LatMinStart / 60
    END AS LatShallow,
CASE 
    WHEN H.TransOrientation = 'D' THEN -(H.LongDegEnd + H.LongMinEnd / 60)
    WHEN H.TransOrientation = 'S' THEN -(H.LongDegStart + H.LongMinStart / 60)
    END AS LonShallow,
H.TransectNum AS Transect,
H.TransLen AS Transect_length,
IsNull(FORMAT (D.QTime,'HH:mm'),FORMAT (H.StartTime,'HH:mm')) AS Time,
H.TimeType AS Time_type,
D.QuadratNum AS Quadrat,
D.GageDepth*0.3048 AS Depth_gauge_m,
D.ChartDepth AS CorDepthM,
D.SubStrate1 AS Substrate1,
D.SubStrate2 AS Substrate2,
CASE 
  	WHEN ',' + A.CanopySpecies1 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.CanopySpecies2 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies1 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies2 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies3 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies4 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies5 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies6 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies7 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies8 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies9 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.UndStySpecies10 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies1 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies2 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies3 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies4 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies5 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + A.TurfSpecies6 + ',' LIKE '%PH%' THEN 1
	ELSE 0 
	END AS PH,
CASE 
	WHEN ',' + A.CanopySpecies1 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.CanopySpecies2 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies1 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies2 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies3 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies4 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies5 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies6 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies7 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies8 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies9 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.UndStySpecies10 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies1 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies2 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies3 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies4 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies5 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + A.TurfSpecies6 + ',' LIKE '%Z%' THEN 1
	ELSE 0 
	END AS ZO
FROM Shellfish_Bio_Urchin.dbo.UrchHeaders H
    INNER JOIN 
    Shellfish_Bio_Urchin.dbo.UrchDensity D ON H.[Key] = D.HKey
    LEFT JOIN
    Shellfish_Bio_Urchin.dbo.UrchHabitat A ON H.[Key] = A.HKey AND D.QuadratNum = A.QuadratNum
WHERE H.LatDegStart IS NOT NULL AND H.LongDegStart IS NOT NULL AND D.ChartDepth IS NOT NULL AND H.Year > 2004 AND H.Species = '6BC' AND H.TransOrientation != 'U' 
ORDER BY H.Year, H.TransectNum, D.QuadratNum;

