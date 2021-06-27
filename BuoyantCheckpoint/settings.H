    dictionary simpleDict = mesh.solutionDict().subDict("SIMPLE");


/*
    if (simpleDict.found("costFunctionPatches"))
    {
        costFunctionPatches.reset
        (
            new wordList(simpleDict.lookup("costFunctionPatches"))
        );
    }else{
        Info << "Warning: Keyword costFunctionPachtes not found in fvSolution/SIMPLE" << endl;
        Info << "Example:\nSIMPLE{\n  costFunctionPatches (inlet outlet);\n}" << endl;
    }
*/
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

    nOptSteps = simpleDict.lookupOrDefault<label>("nOptSteps",5000);
    if (!simpleDict.found("nOptSteps"))
    {
	Info << "Warning: Keyword nOptSteps not found in fvSolution/SIMPLE. Default set to 100" << endl;
    }

    penalty = simpleDict.lookupOrDefault<scalar>("penalty",10);
    if (!simpleDict.found("penalty"))
    {
	Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 10" << endl;
    }
	    
    porosity_f = simpleDict.lookupOrDefault<scalar>("porosity_f",0.0);
    if (!simpleDict.found("porosity_f"))
    {
	Info << "Warning: Keyword porosity_f not found in fvSolution/SIMPLE. Default set to 0.0" << endl;
    }
	    
    porosity_s = simpleDict.lookupOrDefault<scalar>("porosity_s",10e5);
    if (!simpleDict.found("porosity_s"))
    {
	Info << "Warning: Keyword porosity_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
    }
    
    maxPiggyLoop = simpleDict.lookupOrDefault<label>("maxPiggyLoop",100);
    if (!simpleDict.found("maxPiggyLoop"))
    {
	Info << "Warning: Keyword maxPiggyLoop not found in fvSolution/SIMPLE. Default set to 100" << endl;
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
    
    maxMMAiter = simpleDict.lookupOrDefault<label>("maxMMAiter",8);
    if (!simpleDict.found("maxMMAiter"))
    {
	Info << "Warning: Keyword maxMMAiter not found in fvSolution/SIMPLE. Default set to 8" << endl;
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
    
    if (MMALoop == 0) {
	    forAll(designSpaceCells,i){
		const label j = designSpaceCells[i];
		dfdeta[j] = -mesh.V()[j]*designSize/designVolume;
	    }
	    dfdeta.write();
		Info << "Sum dfdeta = " << gSum(dfdeta) << endl;
		
	    
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

