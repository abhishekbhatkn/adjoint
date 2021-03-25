    dictionary simpleDict = mesh.solutionDict().subDict("SIMPLE");

    autoPtr<wordList> costFunctionPatches;
    costFunctionPatches.reset(new wordList(0));

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
	    
    alpha_f = simpleDict.lookupOrDefault<scalar>("alpha_f",0.1);
    if (!simpleDict.found("alpha_f"))
    {
	Info << "Warning: Keyword alpha_f not found in fvSolution/SIMPLE. Default set to 0.1" << endl;
    }
	    
    alpha_s = simpleDict.lookupOrDefault<scalar>("alpha_s",10e5);
    if (!simpleDict.found("alpha_s"))
    {
	Info << "Warning: Keyword alpha_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
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
    
    maxMMA = simpleDict.lookupOrDefault<label>("maxMMA",15);
    if (!simpleDict.found("maxMMA"))
    {
	Info << "Warning: Keyword maxMMA not found in fvSolution/SIMPLE. Default set to 15" << endl;
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
    
    IOdictionary transportProperties
    (
        IOobject
        (
	    "transportProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
        )
    );

    DTs = dimensionedScalar(transportProperties.lookup("DTs"));

    DTf = dimensionedScalar(transportProperties.lookup("DTf"));
    
    n = 0.0;
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	n += mesh.V()[j];
    }
    Foam::reduce(n,sumOp<scalar>());
    
    forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	dfdeta[j] = mesh.V()[j]/n;
    }
    dfdeta.write();
        Info << "Sum dfdeta = " << gSum(dfdeta) << endl;
