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

    penalty = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("penalty",0.6);
    if (!mesh.solutionDict().subDict("SIMPLE").found("penalty"))
	    {
		Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 0.6" << endl;
	    }
	    
    alpha_f = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_f",0.1);
    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_f"))
	    {
		Info << "Warning: Keyword alpha_f not found in fvSolution/SIMPLE. Default set to 0.1" << endl;
	    }
	    
    alpha_s = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_s",10e5);
    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_s"))
	    {
		Info << "Warning: Keyword alpha_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
	    }
	        
    volumeConstraint = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("volumeConstraint",0.4);
    if (!mesh.solutionDict().subDict("SIMPLE").found("volumeConstraint"))
    {
	Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
    }
   
    nOptSteps = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<label>("nOptSteps",100);
    if (!mesh.solutionDict().subDict("SIMPLE").found("nOptSteps"))
    {
	Info << "Warning: Keyword nOptSteps not found in fvSolution/SIMPLE. Default set to 100" << endl;
    }
    
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
