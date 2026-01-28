SELECT 
B.HKey AS HKey,
B.BeltTransect,
B.GaugeDepth,
B.CorDepthM,
B.Substrate1,
B.Sub1Pct,
B.Substrate2,
B.Sub2Pct,
B.Substrate3,
B.Sub3Pct,
B.PA_RSU AS RSU,
B.PA_GSU AS GSU,
B.PA_PSU AS PSU,
B.PA_RSC AS RSC,
B.PA_PT AS PT,
B.PA_PH AS PH,
B.PA_ZM AS ZM
FROM [2 - BeltTransect] B;


