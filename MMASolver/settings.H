    dictionary simpleDict = mesh.solutionDict().subDict("SIMPLE");
  
    const Foam::wordList designSpaceZones(
            simpleDict.lookupOrDefault<Foam::wordList>("designSpace",Foam::wordList())
    );

    if(designSpaceZones.size()>0){
        forAll(designSpaceZones, i){
            const label cellZoneID = mesh.cellZones().findZoneID(designSpaceZones[i]);
            designSpaceCells.append( mesh.cellZones()[cellZoneID] );
        }
    }else{ // add all cells
        forAll(eta, i){
            designSpaceCells.append( i );
        }
    }

    optEpsilon = simpleDict.lookupOrDefault<scalar>("optTolerance",1e-3);
    if (!simpleDict.found("optTolerance"))
    {
	Info << "Warning: Keyword optTolerance not found in fvSolution/SIMPLE. Default set to 1e-3" << endl;
    }
    
    volumeConstraint = simpleDict.lookupOrDefault<scalar>("volumeConstraint",0.4);
    if (!simpleDict.found("volumeConstraint"))
    {
	Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
    }

    factorT = simpleDict.lookupOrDefault<scalar>("factorT",0.0);
    if (!simpleDict.found("factorT"))
    {
	Info << "Warning: Keyword factorT not found in fvSolution/SIMPLE. Default set to 0.0" << endl;
    }
    factorP = simpleDict.lookupOrDefault<scalar>("factorP",1.0);
    if (!simpleDict.found("factorP"))
    {
	Info << "Warning: Keyword factorP not found in fvSolution/SIMPLE. Default set to 1.0" << endl;
    }
    
    minTempCost = simpleDict.lookupOrDefault<scalar>("minTempCost",0.0);
    if (!simpleDict.found("minTempCost"))
    {
	Info << "Warning: Keyword minTempCost not found in fvSolution/SIMPLE. Default set to 0.0" << endl;
    }
    
    maxTempCost = simpleDict.lookupOrDefault<scalar>("maxTempCost",1.0);
    if (!simpleDict.found("maxTempCost"))
    {
	Info << "Warning: Keyword maxTempCost not found in fvSolution/SIMPLE. Default set to 1.0" << endl;
    }
    
    minPrCost = simpleDict.lookupOrDefault<scalar>("minPrCost",0.0);
    if (!simpleDict.found("minPrCost"))
    {
	Info << "Warning: Keyword minPrCost not found in fvSolution/SIMPLE. Default set to 0.0" << endl;
    }
    
    maxPrCost = simpleDict.lookupOrDefault<scalar>("maxPrCost",1.0);
    if (!simpleDict.found("maxPrCost"))
    {
	Info << "Warning: Keyword maxPrCost not found in fvSolution/SIMPLE. Default set to 1.0" << endl;
    }
    
    MMALoop = simpleDict.lookupOrDefault<label>("MMALoop",0);
    if (!simpleDict.found("MMALoop"))
    {
	Info << "Warning: Keyword MMALoop not found in fvSolution/SIMPLE. Default set to 0" << endl;
    }   
    
    asyminit = simpleDict.lookupOrDefault<scalar>("asyminit",0.5);
    if (!simpleDict.found("asyminit"))
    {
	Info << "Warning: Keyword asyminit not found in fvSolution/SIMPLE. Default set to 0.5" << endl;
    }

    asymdec = simpleDict.lookupOrDefault<scalar>("asymdec",0.7);
    if (!simpleDict.found("asymdec"))
    {
	Info << "Warning: Keyword asymdec not found in fvSolution/SIMPLE. Default set to 0.7" << endl;
    }
    
    asyminc = simpleDict.lookupOrDefault<scalar>("asyminc",1.4);
    if (!simpleDict.found("asyminc"))
    {
	Info << "Warning: Keyword asyminc not found in fvSolution/SIMPLE. Default set to 1.4" << endl;
    }
     
    designVolume = 0.0;
    designSize = designSpaceCells.size();
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	designVolume += mesh.V()[j];
    }
    Foam::reduce(designVolume,sumOp<scalar>());
    Foam::reduce(designSize,sumOp<label>());
    
    
    
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	dfdeta[j] = -mesh.V()[j]*designSize/designVolume;
    }
    dfdeta.write();
    Info << "Sum dfdeta = " << gSum(dfdeta) << endl;
    
    if (MMALoop == 0) {	
	    const Foam::wordList solidSpaceZones(
		    mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<Foam::wordList>("solidSpace",Foam::wordList())
	    );
	    forAll(solidSpaceZones, i){
		    const label cellZoneID = mesh.cellZones().findZoneID(solidSpaceZones[i]);
		    solidSpaceCells.append( mesh.cellZones()[cellZoneID] );
		}

	    forAll(solidSpaceCells,i){
	    	const label j = solidSpaceCells[i];
		eta[j] = 0.0;
	    }
	    eta.write();
    }
