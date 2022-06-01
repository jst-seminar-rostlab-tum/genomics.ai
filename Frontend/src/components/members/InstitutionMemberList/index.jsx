import React from 'react';
import LoadingList from 'components/general/LoadingList';
import InstitutionMemberCard from 'components/members/InstitutionMemberCard';

/**
 * @param {isLoading} whether the list is loading
 * @param {members} the members to display
 * @param {nextToNameBuilder} will be called with the member
 *                            and should return a component to display next to the name
 * @param {trailingBuilder} will be called with the member
 *                         and should return a component to display at the end of the card
 * @returns list of members (or loading indicator)
 */
function InstitutionMemberList({
  isLoading, members, nextToNameBuilder, trailingBuilder,
}) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={members}
      cardBuilder={(member) => (
        <InstitutionMemberCard
          member={member}
          nextToNameBuilder={nextToNameBuilder}
          trailingBuilder={trailingBuilder}
        />
      )}
      noElementsMessage="No members."
    />
  );
}

export default InstitutionMemberList;
