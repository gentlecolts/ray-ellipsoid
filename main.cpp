#include <cstdlib>
#include <iostream>
#include <SDL/SDL.h>
#include <cmath>
#include <algorithm>
#include <chrono>
using namespace std;

const float PI=4*atan(1.0);
auto start=
	chrono::system_clock::now().time_since_epoch() /
	chrono::milliseconds(1);

struct vec3d{
	float x,y,z;

	vec3d():vec3d(0,0,0){}
	vec3d(float x0,float y0,float z0):x(x0),y(y0),z(z0){
		//nothing to do here
	}
};
#define dot(v1,v2) (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)

struct ray{
	vec3d ori,dir;
	ray():ray(0,0,0,0,0,0){}
	ray(vec3d& v,vec3d& d):ray(v.x,v.y,v.z,d.x,d.y,d.z){}
	ray(float x0,float y0,float z0,float dx,float dy,float dz):ori(x0,y0,z0),dir(dx,dy,dz){
	}
};

class ellipsoid{
public:
	vec3d pos;
	float x,y,z;

	ellipsoid():pos(0,0,0),x(1),y(1),z(1){}

	bool intersects(ray& r,float &t0,float& t1){
		/*
		P=a point on the ellipsoid
			[	1/a		0		0	]
		M=	[	0		1/b		0	]
			[	0		0		1/c	]
		such that the ellipsoid is defined by
		((x-c1)/a)^2 + ((y-c2)/b)^2 + ((z-c3)/c)^2 = 1

		thus this problem may be expressed via:
		|(P-C)*M|^2=1

		since (A dot A) = |A|^2 we can distribute M, then substitute
		(P*M-C*M) dot (P*M-C*M) = 1

		define the following:
		P=P0+t*v
		P1=P0*M-C*M=(P0-C)*M
		v1=v*M

		substituting gives
		(P0*M + t*v*M - C*M) dot (P0*M + t*v*M - C*M) = 1
		(P0*M - C*M + t*v*M) dot (P0*M - C*M + t*v*M) = 1
		(P1 + t*v1) dot (P1 + t*v1) = 1

		finally, we get
		|P1|^2 + 2*t*(P1 dot v1) + t^2*|v1|^2 = 1

		which may be solved as a quadratic where
		a=|v1|^2
		b=2*(P1 dot v1)
		c=|P1|^2 - 1
		*/
		const vec3d P1(
			(r.ori.x-pos.x)/x,
			(r.ori.y-pos.y)/y,
			(r.ori.z-pos.z)/z
		);
		const vec3d v1(
			r.dir.x/x,
			r.dir.y/y,
			r.dir.z/z
		);

		const float
			a=dot(v1,v1),
			b=2*dot(P1,v1),
			c=dot(P1,P1)-1;

		const float
			eq=b*b-4*a*c,
			rt=sqrt(eq);
		/*
		const float
			r1=(-b-rt)/(2*a),
			r2=(-b+rt)/(2*a);
		t0=min(r1,r2);
		t1=max(r1,r2);
		/*/
		//since in this context, (a) must always be positive, t0 will always be the lower value
		t0=(-b-rt)/(2*a);
		t1=(-b+rt)/(2*a);
		//*/

		return (eq>=0) & ((t0>=0) | (t1>0));
	}
} elip;

unsigned int counter=0;

const float defCamDist=-3;
struct camera{
	vec3d pos,right,up,forward;

	camera():pos(0,0,defCamDist),right(1,0,0),up(0,1,0),forward(0,0,1){}
	camera(SDL_Surface* screen,float fov):
		pos(0,0,defCamDist),right(1,0,0),
		up(0,right.x*float(screen->h)/screen->w,0),
		forward(0,0,right.x/tan(fov/2))//tan(fov/2)=x/z
	{}

	ray getRay(float perX,float perY){
		/*
		ray goes from
			camera origin
		in the direction
			forward + (2*perX-1)*right + (2*perY-1)*up
		*/

		float
			px=2*perX-1,
			py=2*perY-1,
			dirx=forward.x+px*right.x+py*up.x,
			diry=forward.y+px*right.y+py*up.y,
			dirz=forward.z+px*right.z+py*up.z,
			norm=1/sqrt(dirx*dirx+diry*diry+dirz*dirz);
		dirx*=norm;
		diry*=norm;
		dirz*=norm;

		if(!(counter=(counter+1)%1000)){
			//printf("casting ray from (%f,%f,%f) to (%f,%f,%f)\n",pos.x,pos.y,pos.z,dirx,diry,dirz);
		}
		return ray(pos.x,pos.y,pos.z,dirx,diry,dirz);
	}
};

void render(SDL_Surface* screen){
	const int w=screen->w,h=screen->h;
	uint32_t* pix=(uint32_t*)(screen->pixels);

	uint32_t defCol=0x00bfbfbf;

	camera cam(screen,90);

	auto timenow=
		std::chrono::system_clock::now().time_since_epoch() /
		std::chrono::milliseconds(1);
	double
		seconds_since_start = double(timenow-start)/1000,
		angle=2*PI*seconds_since_start;

	//cout<<seconds_since_start<<endl;
	//cam.pos.z=-20+15*sin(angle);

	elip.x=2+1*sin(angle);
	elip.y=1+0.25*cos(angle);
	elip.z=1.5+cos(angle);

	float t0,t1;
	for(int i=0;i<w*h;i++){
		ray camray=cam.getRay(float(i%w)/w,float(i/w)/h);
		bool intr=elip.intersects(camray,t0,t1);
		if(!counter){
			//printf("intersection from: %f to %f\n",t0,t1);
		}

		vec3d negCamRay(
			-camray.dir.x,
			-camray.dir.y,
			-camray.dir.z
		);

		vec3d intersect(
			camray.ori.x+t0*camray.dir.x,
			camray.ori.y+t0*camray.dir.y,
			camray.ori.z+t0*camray.dir.z
		);

		vec3d norm(
			(intersect.x-elip.pos.x)*2/(elip.x*elip.x),
			(intersect.y-elip.pos.y)*2/(elip.y*elip.y),
			(intersect.z-elip.pos.z)*2/(elip.z*elip.z)
		);
		float dirnorm=1/sqrt(
			norm.x*norm.x+
			norm.y*norm.y+
			norm.z*norm.z
		);
		norm.x*=dirnorm;
		norm.y*=dirnorm;
		norm.z*=dirnorm;

		//float p=abs(dot(norm,negNormDir));
		float p=1-acos(dot(norm,negCamRay))/PI;

		/*
		uint32_t
			r=0xe5*p,
			g=0x5c*p,
			b=0x69*p,
			color=(r<<16)|(g<<8)|b;

		pix[i]=intr?color:defCol;
		/*/
		uint32_t
			r=0xe5*p+(1-p)*((defCol>>16)&0xff),
			g=0x5c*p+(1-p)*((defCol>>8)&0xff),
			b=0x69*p+(1-p)*(defCol&0xff),
			color=(r<<16)|(g<<8)|b;

		pix[i]=intr?color:defCol;
		//*/
	}
}

int main ( int argc, char** argv ){
	// initialize SDL video
	if ( SDL_Init( SDL_INIT_VIDEO ) < 0 ){
		printf( "Unable to init SDL: %s\n", SDL_GetError() );
		return 1;
	}

	// make sure SDL cleans up before exit
	atexit(SDL_Quit);

	// create a new window
	SDL_Surface* screen = SDL_SetVideoMode(640, 480, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
	//SDL_Surface* screen = SDL_SetVideoMode(1280, 720, 32, SDL_HWSURFACE|SDL_DOUBLEBUF);
	if ( !screen ){
		printf("Unable to set video: %s\n", SDL_GetError());
		return 1;
	}

	// program main loop
	bool done = false;
	while (!done){
		// message processing loop
		SDL_Event event;
		while (SDL_PollEvent(&event)){
			switch (event.type){
			case SDL_QUIT:
				done = true;
				break;
			case SDL_KEYDOWN:
				// exit if ESCAPE is pressed
				if (event.key.keysym.sym == SDLK_ESCAPE)
					done = true;
				break;
			} // end switch
		} // end of message processing

		// DRAWING STARTS HERE
		render(screen);
		// DRAWING ENDS HERE

		// finally, update the screen :)
		SDL_Flip(screen);
	} // end main loop
}
