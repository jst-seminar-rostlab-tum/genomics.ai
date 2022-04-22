import React from 'react';
import LoadingList from 'components/general/LoadingList';
import MemberCard from 'components/members/MemberCard';

function MemberList({ isLoading, members }) {
  return (
    <LoadingList
      isLoading={isLoading}
      elements={members}
      cardBuilder={(member) => <MemberCard member={member} />}
      noElementsMessage="No members."
    />
  );
}

export default MemberList;
