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

    List<label> designSpaceCells;
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

    	    scalar volumeConstraint = simpleDict.lookupOrDefault<scalar>("volumeConstraint",0.4);
    	    if (!simpleDict.found("volumeConstraint"))
    	    {
        	Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
    	    }
	    scalar penalty = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("penalty",0.6);
	    if (!mesh.solutionDict().subDict("SIMPLE").found("penalty"))
	    {
		Info << "Warning: Keyword penalty not found in fvSolution/SIMPLE. Default set to 0.6" << endl;
	    }
	    dimensionedScalar alpha_f =dimensionedScalar("alpha_f", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_f",0.0));
	    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_f"))
	    {
		Info << "Warning: Keyword alpha_f not found in fvSolution/SIMPLE. Default set to 0" << endl;
	    }
	    dimensionedScalar alpha_s = dimensionedScalar("alpha_s", dimless/dimTime, mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<scalar>("alpha_s",10e5));
	    if (!mesh.solutionDict().subDict("SIMPLE").found("alpha_s"))
	    {
		Info << "Warning: Keyword alpha_s not found in fvSolution/SIMPLE. Default set to 10e5" << endl;
	    }

