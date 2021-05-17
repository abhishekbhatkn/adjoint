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

    nOptSteps = simpleDict.lookupOrDefault<label>("nOptSteps",100);
    if (!simpleDict.found("nOptSteps"))
    {
	Info << "Warning: Keyword nOptSteps not found in fvSolution/SIMPLE. Default set to 100" << endl;
    }

    penalty = simpleDict.lookupOrDefault<scalar>("penalty",0.6);
    if (!simpleDict.found("penalty"))
    {
	Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 0.6" << endl;
    }
	    
    porosity_f = simpleDict.lookupOrDefault<scalar>("porosity_f",0.1);
    if (!simpleDict.found("porosity_f"))
    {
	Info << "Warning: Keyword porosity_f not found in fvSolution/SIMPLE. Default set to 0.1" << endl;
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
    
    optEpsilon = simpleDict.lookupOrDefault<scalar>("optTolerance",5e-2);
    if (!simpleDict.found("optTolerance"))
    {
	Info << "Warning: Keyword optTolerance not found in fvSolution/SIMPLE. Default set to 5e-2" << endl;
    }
    
    initLoop = simpleDict.lookupOrDefault<label>("initLoop",5);
    if (!simpleDict.found("initLoop"))
    {
	Info << "Warning: Keyword initLoop not found in fvSolution/SIMPLE. Default set to 5" << endl;
    }
    
    volumeConstraint = simpleDict.lookupOrDefault<scalar>("volumeConstraint",0.4);
    if (!simpleDict.found("volumeConstraint"))
    {
	Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
    }
    
    maxoutit = simpleDict.lookupOrDefault<label>("maxoutit",8);
    if (!simpleDict.found("maxoutit"))
    {
	Info << "Warning: Keyword maxoutit not found in fvSolution/SIMPLE. Default set to 8" << endl;
    }
    	
    maxInnerLoop = simpleDict.lookupOrDefault<label>("maxInnerLoop",15);
    if (!simpleDict.found("maxInnerLoop"))
    {
	Info << "Warning: Keyword maxLoop not found in fvSolution/SIMPLE. Default set to 15" << endl;
    }
    
    factorT = simpleDict.lookupOrDefault<scalar>("factorT",0.0);
    factorP = simpleDict.lookupOrDefault<scalar>("factorP",1.0);
    
    // reference rho
	refRho = 1027;

    // reference Cp
	Cp = 4195;
    
    designVolume = 0.0;
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	designVolume += mesh.V()[j];
    }
    Foam::reduce(designVolume,sumOp<scalar>());
    
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	dfdeta[j] = mesh.V()[j]/designVolume;
    }
    dfdeta.write();
        Info << "Sum dfdeta = " << gSum(dfdeta) << endl;
