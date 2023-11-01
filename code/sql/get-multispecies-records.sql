SELECT 
	H.TRIP_ID AS Trip_id,
	S.FE_PARENT_EVENT_ID AS HKey,
	H.Survey AS Survey,
	H.Year AS Year,
	H.Month AS Month,
	H.Day AS Day,
	H.StatArea AS Major_stat_area_code,
	H.SubArea AS Stat_subarea_code,
	H.latStart AS LatDeep,
	H.lonStart AS LonDeep,
	H.latEnd AS LatShallow,
	H.lonEnd AS LonShallow,
	H.TransectSite AS Transect,
	IsNull(FORMAT (D.Time,'HH:mm'),FORMAT (H.TimeIn,'HH:mm')) AS Time,
	D.Quadrat AS Quadrat,
	D.GaugeDepth*0.3048 AS Depth_gauge_m,
	D.CorDepthM,
	H.QuadratSkipping+1 AS Quadrat_distance,
	CASE 
		WHEN H.QuadratSkipping = 0 THEN 25
		WHEN H.QuadratSkipping = 1 THEN 50
		WHEN H.QuadratSkipping = 2 THEN 75
		WHEN H.QuadratSkipping = 3 THEN 100
		WHEN H.QuadratSkipping = 4 THEN 125
		END AS Transect_length,
	cast(D.SubCode1 AS INTEGER) AS Substrate1,
	cast(D.SubCode2 AS INTEGER) AS Substrate2,
	CASE WHEN ',' + Q.ALGAE + ',' LIKE '%, PH,%' THEN 1 ELSE 0 END AS PH,
	CASE WHEN ',' + Q.ALGAE + ',' LIKE '%, ZO,%' THEN 1 ELSE 0 END AS ZO
FROM SFBioSQL.dbo.vw_ms1_headers H
	INNER JOIN 
	SFBioSQL.dbo.vw_ms2_density D ON H.TRIP_ID = D.TRIP_ID AND H.TransectSite = D.TransectSite
	INNER JOIN 
	SFBioSQL.dbo.vw_ms2b_summarycounts C ON H.TRIP_ID = C.TRIP_ID AND H.TransectSite = C.TransectSite AND D.Quadrat = C.Quadrat
	INNER JOIN
	SFBioSQL.dbo.vw_ms3_algae2 Q ON H.TRIP_ID = Q.TRIP_ID AND H.TransectSite = Q.TransectSite AND D.Quadrat = Q.Quadrat
	INNER JOIN
	SFBioSQL.dbo.vw_ms2a_subcat S ON H.TRIP_ID = S.TRIP_ID AND H.TransectSite = S.TransectSite AND D.Quadrat = S.Quadrat
WHERE H.latStart IS NOT NULL AND H.lonStart IS NOT NULL AND D.CorDepthM IS NOT NULL 
ORDER BY H.Year, H.TransectSite, D.Quadrat
;
