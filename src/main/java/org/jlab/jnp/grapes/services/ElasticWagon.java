package org.jlab.jnp.grapes.services;

import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.LorentzVector;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;

/**
 * 
 * DVCS Skimming
 *
 * @author fxgirod
 */

public class ElasticWagon extends BeamTargetWagon {

	public ElasticWagon() {
		super("ElasticWagon","fxgirod","0.0");
	}

	public double Vangle(Vector3 v1, Vector3 v2){
		double res=0;
		double l1 = v1.mag();
		double l2 = v2.mag();
		if( l1*l2 > 0)res = Math.toDegrees( Math.acos( v1.dot(v2)/(l1*l2) ) );
		return res;
	}

	@Override
	public boolean processDataEvent(Event event, SchemaFactory factory) {
		LorentzVector VB = new LorentzVector(0,0,beamEnergy,beamEnergy);
		LorentzVector VT = new LorentzVector(0,0,0,targetMass);

		Bank RecPart = new Bank(factory.getSchema("REC::Particle"));
		event.read(RecPart);

		boolean hasElastic = false;
		if(RecPart.getRows()>1 ){
			int n_e = 0;
			int n_p = 0;
			for (int ii = 0; ii < RecPart.getRows() ; ii++) {
				int is_pid = RecPart.getInt("pid", ii);
				int stat   = Math.abs(RecPart.getShort("status", ii));
				if(stat>2000  && stat<4000 && is_pid==11)n_e++;
				if(stat>2000  && stat!=4000 && is_pid==2212)n_p++;
			}
			boolean is_candidate = n_e*n_p>0;
			if(is_candidate){
				int[] e_ind   = new int[n_e];
				int[] p_ind   = new int[n_p];
				n_e = 0;
				n_p = 0;
				for (int ii = 0; ii < RecPart.getRows() ; ii++) {
					int is_pid = RecPart.getInt("pid", ii);
					int stat   = Math.abs(RecPart.getShort("status", ii));
					if(stat>2000  && stat<4000 && is_pid==11){e_ind[n_e]=ii;n_e++;}
					if(stat>2000  && stat!=4000 && is_pid==2212){p_ind[n_p]=ii;n_p++;}
				}
				for (int ie = 0; ie < n_e && !hasElastic; ie++) {
					double e_px  = RecPart.getFloat("px", e_ind[ie]);
					double e_py  = RecPart.getFloat("py", e_ind[ie]);
					double e_pz  = RecPart.getFloat("pz", e_ind[ie]);

					double e_mom = Math.sqrt(e_px*e_px+e_py*e_py+e_pz*e_pz);
					if( e_mom>0.1*beamEnergy ){
						LorentzVector VE = new LorentzVector(e_px,e_py,e_pz,e_mom);
						for (int ip = 0; ip < n_p && !hasElastic; ip++) {
							double p_px  = RecPart.getFloat("px", p_ind[ip]);
							double p_py  = RecPart.getFloat("py", p_ind[ip]);
							double p_pz  = RecPart.getFloat("pz", p_ind[ip]);

							double p_ene = Math.sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz+targetMass*targetMass);
							if( p_ene>0.94358 ){
								LorentzVector VP = new LorentzVector(p_px,p_py,p_pz,p_ene);
								
								LorentzVector Q = new LorentzVector(0,0,0,0);
								Q.add(VB);
								Q.sub(VE);
								LorentzVector W = new LorentzVector(0,0,0,0);
								W.add(Q);
								W.add(VT);

								double delPhi = Math.toDegrees( VE.phi() - VP.phi() );
								while(delPhi>360)delPhi-=360;
								while(delPhi<  0)delPhi+=360;
								boolean isBack2Back = Math.abs( delPhi-180 ) < 7.5;

								double tanTan = Math.tan(VP.theta())*Math.tan(0.5*VE.theta());
								double rEb = targetMass * (1-tanTan)/tanTan;

								//pure elastic condition
								hasElastic = true
									&& W.mass() < 1.6
									&& isBack2Back
									&& Math.abs( rEb - beamEnergy ) < 0.12*beamEnergy 
									;
								if(false && !hasElastic){
									double radElastElecMom = rEb / ( 1  + rEb * ( 1 - Math.cos(VE.theta() ) ) / targetMass );
									
									LorentzVector Vmiss = new LorentzVector(0,0,0,0);
									Vmiss.add(W);
									Vmiss.sub(VP);
									double missDir = Math.sqrt( Vmiss.px()*Vmiss.px() + Vmiss.py()*Vmiss.py() ) / Vmiss.p();

									hasElastic = true // pick up radiative events
										&& Math.abs( e_mom - radElastElecMom ) < 0.12*radElastElecMom
										&& isBack2Back
										&& missDir < 0.25
										;
								}
							}//minimum proton momentum
						}//loop over proton list
					}//minimum electron momentum
				}//loop over electron list
			}//at least one electron and one proton
		}// at least 2 particles
		return hasElastic;
	}//process event
}//elastic wagon
