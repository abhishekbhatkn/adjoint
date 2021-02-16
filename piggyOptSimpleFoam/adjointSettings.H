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
/*
    scalar volumeConstraint = simpleDict.lookupOrDefault<scalar>("volumeConstraint",0.4);
    if (!simpleDict.found("volumeConstraint"))
    {
        Info << "Warning: Keyword volumeConstraint not found in fvSolution/SIMPLE. Default set to 0.4" << endl;
    }
*/
