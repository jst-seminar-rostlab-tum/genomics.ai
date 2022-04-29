import React, { useState, useEffect } from 'react';
import ListCard from 'components/general/ListCard';
import { CircularProgress } from '@mui/material';
import getMember from 'shared/services/mock/members';

/**
 * @param {int} memberId
 * @param {Function(Object)} nextToNameBuilder will be called with the member
 *                            and should return a component to display next to the name
 * @param {Function(Object)} trailingBuilder will be called with the member
 * @returns the card for the member
 */
function MemberCard({ memberId, nextToNameBuilder, trailingBuilder }) {
  const [isLoading, setIsLoading] = useState(true);
  const [member, setMember] = useState({});

  useEffect(async () => {
    setIsLoading(true);
    setMember(await getMember(memberId));
    setIsLoading(false);
  }, [memberId, setMember, setIsLoading]);

  return (
    <ListCard
      imageURL={member.profilePictureURL}
      enforceImage
      title={isLoading ? '...' : `${member.firstName} ${member.lastName}`}
      description={isLoading ? '...' : member.email}
      nextToTitle={isLoading ? null : nextToNameBuilder(member)}
      trailing={isLoading ? <CircularProgress /> : trailingBuilder(member)}
    />
  );
}

export default MemberCard;
