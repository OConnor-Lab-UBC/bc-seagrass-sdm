SELECT
H.[Key] AS HKey,
H.Source AS Survey,
H.Year AS Year,
H.Month AS Month,
H.Day AS Day,
H.StatArea AS Major_stat_area_code,
H.SubArea AS Stat_subarea_code,
H.LatDegSt + H.LatMinSt / 60 AS LatDeep,
-(H.LongDegSt + H.LongMinSt / 60) AS LonDeep,
H.LatDegEnd + H.LatMinEnd / 60 AS LatShallow,
-(H.LongDegEnd + H.LongMinEnd / 60) AS LonShallow,
H.Transect AS Transect,
H.TransectLength AS Transect_length,
IsNull(FORMAT (D.QTime,'HH:mm'),FORMAT (H.StartTime,'HH:mm')) AS Time,
H.TimeType AS Time_type,
D.QuadratNum AS Quadrat,
D.TransDist AS Transect_dist_from_start,
D.DepthGauge *0.3048 AS Depth_gauge_m,
D.DepthCorM AS CorDepthM,
D.Substrate1,
D.Substrate2,
CASE WHEN ',' + D.Algae + ',' LIKE '%PH%' THEN 1 ELSE 0 END AS PH, 
CASE WHEN ',' + D.Algae + ',' LIKE '%ZO%' THEN 1 ELSE 0 END AS ZO
FROM Shellfish_Bio_Other.dbo.GeoduckHeaders H
    INNER JOIN 
    Shellfish_Bio_Other.dbo.GeoduckDensity D ON H.[Key] = D.HKey
WHERE H.LatDegSt IS NOT NULL AND H.LongDegSt IS NOT NULL AND D.DepthCorM IS NOT NULL AND H.Year > 2004 AND D.QuadratNum > 0
ORDER BY H.Year, H.Transect, D.QuadratNum;
