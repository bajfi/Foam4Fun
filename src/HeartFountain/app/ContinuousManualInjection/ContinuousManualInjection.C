// Copyright (c) 2026 JackLee
//
// This software is released under the MIT License.
// https://opensource.org/licenses/MIT

#include "ContinuousManualInjection.H"
#include "KinematicCloud.H"
#include "basicKinematicCloud.H"
#include "basicKinematicCollidingCloud.H"
#include "bitSet.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

namespace Foam
{
makeInjectionModelType(ContinuousManualInjection, basicKinematicCollidingCloud);

template <>
const word GlobalIOList<Tuple2<vector, vector>>::typeName("Tuple2VectorVector");
} // namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ContinuousManualInjection<CloudType>::ContinuousManualInjection(
  const dictionary& dict, CloudType& owner, const word& modelName
)
  : InjectionModel<CloudType>(dict, owner, modelName, typeName),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positionAxis_(IOobject(
      positionsFile_,
      owner.db().time().constant(),
      owner.mesh(),
      IOobject::MUST_READ,
      IOobject::NO_WRITE,
      IOobject::NO_REGISTER
    )),
    Umag_(Function1<scalar>::New("Umag", this->coeffDict(), &owner.mesh())),
    injectorCells_(positionAxis_.size(), -1),
    injectorTetFaces_(positionAxis_.size(), -1),
    injectorTetPts_(positionAxis_.size(), -1),
    sizeDistribution_(
      distributionModel::New(this->coeffDict().subDict("sizeDistribution"), owner.rndGen())
    ),
    ignoreOutOfBounds_(this->coeffDict().getOrDefault("ignoreOutOfBounds", false)),
    duration_(this->coeffDict().getScalar("duration")),
    parcelsPerSecond_(this->coeffDict().template get<label>("parcelsPerSecond")),
    flowRateProfile_(Function1<scalar>::New("flowRateProfile", this->coeffDict(), &owner.mesh())),
    injectorOrder_(identity(positionAxis_.size()))
{
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    flowRateProfile_->userTimeToTime(time);

    Umag_->userTimeToTime(time);

    updateMesh();

    // Determine volume of particles to inject
    this->volumeTotal_ = flowRateProfile_->integrate(0.0, duration_);
}

template <class CloudType>
Foam::ContinuousManualInjection<CloudType>::ContinuousManualInjection(
  const ContinuousManualInjection<CloudType>& im
)
  : InjectionModel<CloudType>(im),
    positionsFile_(im.positionsFile_),
    positionAxis_(im.positionAxis_),
    Umag_(im.Umag_->clone().ptr()),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    sizeDistribution_(im.sizeDistribution_.clone()),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_.clone().ptr()),
    injectorOrder_(im.injectorOrder_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ContinuousManualInjection<CloudType>::~ContinuousManualInjection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::ContinuousManualInjection<CloudType>::updateMesh()
{
    bitSet keep(positionAxis_.size(), true);
    label nRejected = 0;

    // Set/cache the injector cells
    forAll(positionAxis_, i)
    {
        if (!this->findCellAtPosition(
              injectorCells_[i],
              injectorTetFaces_[i],
              injectorTetPts_[i],
              positionAxis_[i].first(),
              !ignoreOutOfBounds_
            ))
        {
            keep.unset(i);
            ++nRejected;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);
        inplaceSubset(keep, positionAxis_);
    }
}

template <class CloudType>
Foam::scalar Foam::ContinuousManualInjection<CloudType>::timeEnd() const
{
    // Injection is instantaneous - but allow for a finite interval to
    // avoid numerical issues when interval is zero
    return this->SOI_ + duration_;
}

template <class CloudType>
Foam::label Foam::ContinuousManualInjection<CloudType>::parcelsToInject(
  const scalar time0, const scalar time1
)
{
    if (positionAxis_.empty())
    {
        return 0;
    }

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        const scalar t1 = min(time1, duration_);
        const scalar nParcels = max(0.0, (t1 - time0) * parcelsPerSecond_);
        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        Random& rnd = this->owner().rndGen();
        scalar rndPos = rnd.globalPosition(scalar(0), scalar(1));
        if (nParcels - scalar(nParcelsToInject) > rndPos)
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }

    return 0;
}

template <class CloudType>
Foam::scalar Foam::ContinuousManualInjection<CloudType>::volumeToInject(
  const scalar time0, const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_->integrate(time0, min(time1, duration_));
    }

    return 0.0;
}

template <class CloudType>
void Foam::ContinuousManualInjection<CloudType>::setPositionAndCell(
  const label parcelI,
  const label nParcels,
  const scalar time,
  vector& position,
  label& cellOwner,
  label& tetFacei,
  label& tetPti
)
{
    // Random choose a position from the list of positions
    Random& rnd = this->owner().rndGen();
    rnd.shuffle(injectorOrder_);
    const label i = injectorOrder_[parcelI % positionAxis_.size()];

    position = positionAxis_[i].first();
    cellOwner = injectorCells_[i];
    tetFacei = injectorTetFaces_[i];
    tetPti = injectorTetPts_[i];
}

template <class CloudType>
void Foam::ContinuousManualInjection<CloudType>::setProperties(
  const label parcelI,
  const label nParcels,
  const scalar time,
  typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() =
      positionAxis_[injectorOrder_[parcelI % positionAxis_.size()]].second() * Umag_->value(time);

    // set particle diameter
    parcel.d() = sizeDistribution_->sample();
}

template <class CloudType>
bool Foam::ContinuousManualInjection<CloudType>::fullyDescribed() const
{
    return false;
}

template <class CloudType>
bool Foam::ContinuousManualInjection<CloudType>::validInjection(const label)
{
    return true;
}

// ************************************************************************* //
