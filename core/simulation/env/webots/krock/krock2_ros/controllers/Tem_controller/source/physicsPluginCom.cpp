// ========================================================================================================================================================================================

// File:          physicsPluginCom.cpp
// Date:          25-07-2014
// Description:   functions to communicate between controller and physics plugin
// Author:        Robin Thandiackal
// Modifications: -

// ========================================================================================================================================================================================

#include "physicsPluginCom.hpp"

connection::connection()
{

}

bool
connection::communicate2physics(webots::Emitter *emitter, webots::Receiver *receiver, c2p_packet e_packet, p2c_packet &r_packet)
{

	bool a=0,b=0;
	a=send(emitter,e_packet);
	b=receive(receiver,r_packet);
	
	
	/*b=1;
	p2c_packet packet;
	for(int i=0;i<N_SERVOS;i++)
	{
		packet.torques[i] = 0.0;
	}	
	receive(receiver, packet);
	cout << packet.torques[12] << endl;*/


	if(a && b)
	{
		//cerr << "a: " << a << " b: " << b << endl;
		//cout << sizeof(r_packet) << endl;
		//cout << "hallo" << endl;
		//cout << r_packet->torques[12] << endl;
		return 1;
	}	
	else
	{
		return 0;
	}
}

bool
connection::send(webots::Emitter *emitter, c2p_packet packet)
{
	if(sizeof(packet)>=1024)
	{
    		std::cerr << " size of packet: " << sizeof(packet) << std::endl;
		std::cerr << " PACKET SIZE EXCEEDS BUFFER SIZE: LIMITED TO 1024 !!!! " << std::endl;
		return 0;
    	}
    	else
    	{
      		emitter->send(&packet, sizeof(packet));
		return 1;
    	}
}
	
bool
connection::receive(webots::Receiver *receiver, p2c_packet &packet)
{
	p2c_packet *p;

	if (receiver->getQueueLength()>0)
    	{
      		p = (p2c_packet*)(receiver->getData());
		packet = *p;
      		receiver->nextPacket(); // very important!!!!! (deletes the head packet, the next packet becomes new head!)
		//cout << packet->torques[12] << endl;
		return 1;
    	}
	else
	{		
		std::cerr << " NO PACKET RECEIVED FROM PLUGIN" << std::endl;
		return 0;
	}
}


connection::~connection()
{

}
