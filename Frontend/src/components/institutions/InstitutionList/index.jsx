import React from 'react';
import LoadingList from 'components/general/LoadingList';
import InstitutionCard from 'components/institutions/InstitutionCard';

function InstitutionList({
  isLoading, institutions, disableLinks, replaceTrailing,
}) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={institutions}
      cardBuilder={(institution) => (
        <InstitutionCard
          institution={institution}
          replaceTrailing={replaceTrailing}
          disableLink={disableLinks}
        />
      )}
      noElementsMessage="No institutions."
    />
  );
}

export default InstitutionList;
