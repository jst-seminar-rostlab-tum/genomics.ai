import React from 'react';
import LoadingList from 'components/general/LoadingList';
import InstitutionCard from 'components/institutions/InstitutionCard';

function InstitutionList({ isLoading, institutions }) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={institutions}
      cardBuilder={(institution) => <InstitutionCard institution={institution} />}
      noElementsMessage="No institutions."
    />
  );
}

export default InstitutionList;
