SELECT
H.[Key] AS HKey,
H.Source AS Survey,
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
H.TransLen AS Transect_length,
IsNull(FORMAT (D.QTime,'HH:mm'),FORMAT (H.StartTime,'HH:mm')) AS Time,
H.TimeType AS Time_type,
D.QuadratNum AS Quadrat,
CASE 
    WHEN H.[Key] IN (11423, 9014, 9005, 10882, 10881, 8705, 8699, 8698, 8697, 8695, 8691, 8610, 8600, 8592) THEN 1
	WHEN H.SurveyPattern IS NULL THEN 1
	WHEN H.SurveyPattern IN ('MCMC', 'MMMM', 'MMM', 'MCMCM', 'MMMM…MCMC', 'MMMM…CCCC', 'MCCC') THEN 1
	WHEN H.SurveyPattern IN ('MSMS', 'MSCSM', 'MSCS', 'MSCSCS', 'MSCSCSCSCS', '*MSCSCSCSCS', '*MSCSCS', 'MSCSCSCS', '*MMMMMMSCSCS') THEN 2
	WHEN H.SurveyPattern IN ('*MSSCSSCSSCSSCSS', 'MSSCSSCSSCSS') THEN 3
	END AS Quadrat_distance,
D.GageDepth*0.3048 AS Depth_gauge_m,
CASE 
    WHEN H.Year < 2010 THEN (D.GageDepth - D.TideHgt)*0.3048
    WHEN H.Year > 2009 THEN D.ChartDepth
    END AS CorDepthM,	
D.SubStrate1 AS Substrate1,
D.PctSub1 AS PctSubstrate1,
D.SubStrate2 AS Substrate2,
D.PctSub2 AS PctSubstrate2,
D.SubStrate3 AS Substrate3,
D.PctSub3 AS PctSubstrate3,
A.UnderstoryPct,
CASE 
    WHEN ',' + D.Algae + ',' LIKE '%PH%' THEN 1 
	WHEN ',' + D.Algae2 + ',' LIKE '%PH%' THEN 1
	WHEN ',' + D.Algae3 + ',' LIKE '%PH%' THEN 1 
	WHEN ',' + D.Algae4 + ',' LIKE '%PH%' THEN 1
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
    WHEN ',' + D.Algae + ',' LIKE '%Z%' THEN 1 
	WHEN ',' + D.Algae2 + ',' LIKE '%Z%' THEN 1
	WHEN ',' + D.Algae3 + ',' LIKE '%Z%' THEN 1 
	WHEN ',' + D.Algae4 + ',' LIKE '%Z%' THEN 1
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
WHERE H.LatDegStart IS NOT NULL AND H.LongDegStart IS NOT NULL AND D.ChartDepth IS NOT NULL AND H.Year BETWEEN 1993 AND 2023 AND H.Species = '6BB' AND H.TransOrientation != 'U' 
ORDER BY H.Year, H.[Key], D.QuadratNum;


