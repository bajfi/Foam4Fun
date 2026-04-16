// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#pragma once

#include "InjectionModel.H"
#include "Switch.H"
#include "distributionModel.H"
#include "GlobalIOList.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ContinuousManualInjection Declaration
\*---------------------------------------------------------------------------*/

template <class CloudType>
class ContinuousManualInjection : public InjectionModel<CloudType>
{
    // Private data

    //- Name of file containing positions data
    const word positionsFile_;

    //- Global list of position and axis for each injector
    GlobalIOList<Tuple2<vector, vector>> positionAxis_;

    //- Magnitude of velocity to inject particles with
    autoPtr<Function1<scalar>> Umag_;

    //- List of cell labels corresponding to injector positions
    labelList injectorCells_;

    //- List of tetFace labels corresponding to injector positions
    labelList injectorTetFaces_;

    //- List of tetPt labels corresponding to injector positions
    labelList injectorTetPts_;

    //- Parcel size distribution model
    const autoPtr<distributionModel> sizeDistribution_;

    //- Flag to suppress errors if particle injection site is out-of-bounds
    Switch ignoreOutOfBounds_;

    //- Duration to injection particles
    scalar duration_;

    //- Number of parcels to introduce per second []
    const label parcelsPerSecond_;

    //- Flow rate profile relative to SOI []
    autoPtr<Function1<scalar>> flowRateProfile_;

    //- Order to process the injectors
    labelList injectorOrder_;

  public:
    //- Runtime type information
    TypeName("continuousManualInjection");

    // Constructors

    //- Construct from dictionary
    ContinuousManualInjection(const dictionary& dict, CloudType& owner, const word& modelName);

    //- Construct copy
    ContinuousManualInjection(const ContinuousManualInjection<CloudType>& im);

    //- Construct and return a clone
    virtual autoPtr<InjectionModel<CloudType>> clone() const
    {
        return autoPtr<InjectionModel<CloudType>>(new ContinuousManualInjection<CloudType>(*this));
    }

    //- Destructor
    virtual ~ContinuousManualInjection();

    // Member Functions

    //- Set injector locations when mesh is updated
    virtual void updateMesh();

    //- Return the end-of-injection time
    scalar timeEnd() const;

    //- Number of parcels to introduce relative to SOI
    virtual label parcelsToInject(const scalar time0, const scalar time1);

    //- Volume of parcels to introduce relative to SOI
    virtual scalar volumeToInject(const scalar time0, const scalar time1);

    // Injection geometry

    //- Set the injection position and owner cell, tetFace and tetPt
    virtual void setPositionAndCell(
      const label parcelI,
      const label nParcels,
      const scalar time,
      vector& position,
      label& cellOwner,
      label& tetFacei,
      label& tetPti
    );

    //- Set the parcel properties
    virtual void setProperties(
      const label parcelI,
      const label nParcels,
      const scalar time,
      typename CloudType::parcelType& parcel
    );

    //- Flag to identify whether model fully describes the parcel
    virtual bool fullyDescribed() const;

    //- Return flag to identify whether or not injection of parcelI is
    //  permitted
    virtual bool validInjection(const label parcelI);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
