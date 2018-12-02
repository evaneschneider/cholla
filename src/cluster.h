/*! \file cluster.h
 *  \brief Declarations of the cluster class. */

#ifndef CLUSTER_H
#define CLUSTER_H

#include"global.h"

extern float S99_table[1000][3];
extern float cluster_list[N_CL][5];

/*! \class Cluster
 *  \brief Class to create a supernova cluster. */

class Cluster
{
  public:

    /*! \var id
     *  \brief ID of the cluster */
    int id;

    /*! \var flag_on
     *  \brief Flag to tell whether the cluster is on or off.
               0 is off, 1 is on (first step), 2 is on */
    int flag_on;

    /*! \var mass
     *  \brief Mass of the cluster */
    int mass;

    /*! \var time 
     *  \brief Time the cluster has been on */
    Real time;

    /*! \var SF_cl
     *  \brief Total cumulative star formation as of this point
         in the cluster list */
    Real SF_cl;

    /*! \var x_pos
     *  \brief Global x_position of the cluster */
    Real x_pos;

    /*! \var y_pos
     *  \brief Global y_position of the cluster */
    Real y_pos;

    /*! \var z_pos
     *  \brief Global z_position of the cluster */
    Real z_pos;

    /*! \var r_pos
     *  \brief Global radial position of the cluster */
    Real r_pos;
    
    /*! \var phi_pos
     *  \brief Global phi position of the cluster */
    Real phi_pos;

    /*! \var R
     *  \brief Radius of the cluster (in kpc) */
    Real R_cl;

    /*! \var V_cl
    *  \brief Volume of the cluster */
    Real V_cl;

    /*! \var M_dot
     *  \brief Current M_dot from S99 fluxes */
    Real M_dot;

    /*! \var E_dot
     *  \brief Current E_dot from S99 fluxes */
    Real E_dot;

    /*! \fn Initialize(void)
     *  \brief Set the initial cluster variables */
    void Initialize(void);

    /*! \fn Rotate(Real dt)
     *  \brief Rotate the cluster position given a time step */
    void Rotate(Real dt);

    /*! \fn Switch(Real t)
     *  \brief Turn the clusters on and off according to desired star formation */
    void Switch(Real t);

    void Get_S99_Fluxes(void);
    
};

/* \fn void Load_Cluster_List()
 * \brief Load the table of clusters into memory. */
void Load_Cluster_List();

/* \fn void Load_S99_Table()
 * \brief Load the table of M_dot and E_dot from Starburst 99. */
void Load_S99_Table();


#endif //CLUSTER_H

