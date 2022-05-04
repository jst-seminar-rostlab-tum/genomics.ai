import React from 'react';
import ListCard from 'components/general/ListCard';

/**
 * @param {int} memberId
 * @param {Function(Object)} nextToNameBuilder will be called with the member
 *                            and should return a component to display next to the name
 * @param {Function(Object)} trailingBuilder will be called with the member
 * @returns the card for the member
 */
function MemberCard({
  member, nextToNameBuilder, trailingBuilder, overrideProfilePicture = null,
}) {
  return (
    <ListCard
      imageURL={overrideProfilePicture || member.avatarUrl}
      enforceImage
      title={`${member.firstName} ${member.lastName}`}
      description={member.email}
      nextToTitle={nextToNameBuilder ? nextToNameBuilder(member) : null}
      trailing={(trailingBuilder || (() => null))(member)}
    />
  );
}

export default MemberCard;
