{
	forAll(designSpaceCells,i){
		const label j = designSpaceCells[i];
		scalar valEta = eta[j];
		//scalar val = alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty);
		alpha[j] = (alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty));
	}

	//DT = (DTs+(DTf-DTs)*eta*(1.0+penalty)/(eta + penalty));

}
