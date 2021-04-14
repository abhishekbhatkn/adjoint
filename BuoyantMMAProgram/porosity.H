{
	forAll(designSpaceCells,i){
		const label j = designSpaceCells[i];
		scalar valEta = eta[j];
		//scalar val = alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty);
		porosity[j] = (porosity_s+(porosity_f-porosity_s)*valEta*(1.0+penalty)/(valEta + penalty));
	}
	porosity.correctBoundaryConditions();
	
	//alphaTemp[j] = (porosity_s+(porosity_f-porosity_s)*valEta*(1.0+penalty)/(valEta + penalty));
	//DT = (DTs+(DTf-DTs)*eta*(1.0+penalty)/(eta + penalty));

}
