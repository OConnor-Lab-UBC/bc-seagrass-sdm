SELECT 
H.Survey AS Survey, 
H.HKey,
Q.Quadrat, 
S.SpType,
S.Species AS SpeciesCode
FROM ([1-Headers] H
	INNER JOIN [2-Quadrat] Q ON H.HKey = Q.HKey) 
	INNER JOIN [3-SpeciesObservations] S ON (Q.Quadrat = S.Quadrat) AND (Q.HKey = S.HKey)
ORDER BY H.Year, H.Transect, Q.Quadrat;

