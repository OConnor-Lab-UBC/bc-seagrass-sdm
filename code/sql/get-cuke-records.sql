SELECT
H.[Key] AS HKey,
H.Source AS Survey,
H.Year AS Year,
H.Month AS Month,
H.Day AS Day,
H.StatArea AS Major_stat_area_code,
H.SubArea AS Stat_subarea_code,
H.LatDegDeep + H.LatMinDeep / 60 AS LatDeep,
-(H.LongDegDeep + H.LongMinDeep / 60) AS LonDeep,
H.LatDegShallow + H.LatMinShallow / 60 AS LatShallow,
-(H.LongDegShallow + H.LongMinShallow / 60) AS LonShallow,
H.Transect,
H.TransLen AS Transect_length,
IsNull(FORMAT (D.QTime,'HH:mm'),FORMAT (H.TimeStart,'HH:mm')) AS Time,
H.TimeType AS Time_type,
D.QuadratNum AS Quadrat,
D.GaugeDepthFeet*0.3048 AS Depth_gauge_m,
D.ChartDepth AS CorDepthM,
H.MaxChartDepth,
MaxGageDepth*0.3048 AS Max_gauge_m,
SubstrateRaw1 AS Substrate1,
SubstrateRaw2 AS Substrate2,
CASE WHEN ',' + D.Algae1 + ',' LIKE '%PH%' THEN 1 WHEN ',' + D.Algae2 + ',' LIKE '%PH%' THEN 1 ELSE 0 END AS PH,
CASE WHEN ',' + D.Algae1 + ',' LIKE '%ZO%' THEN 1 WHEN ',' + D.Algae2 + ',' LIKE '%ZO%' THEN 1 ELSE 0 END AS ZO
FROM Shellfish_Bio_Other.dbo.SeaCukeHeaders H
    INNER JOIN 
    Shellfish_Bio_Other.dbo.SeaCukeDensities D ON H.[Key] = D.HKey
WHERE H.LatDegDeep IS NOT NULL AND H.LongDegDeep IS NOT NULL AND D.ChartDepth IS NOT NULL AND H.Year > 1999
ORDER BY H.Year, H.[Key], H.Transect, D.QuadratNum;
