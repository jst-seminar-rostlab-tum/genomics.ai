import React from 'react';
import { useHistory } from 'react-router-dom';
import InstitutionMemberListCard from 'components/general/InstitutionMemberListCard';

/**
 * @param {int} memberId
 * @param {Function(Object)} nextToNameBuilder will be called with the member
 *                            and should return a component to display next to the name
 * @param {Function(Object)} trailingBuilder will be called with the member
 * @returns the card for the member
 */
function InstitutionMemberCard({
  member, nextToNameBuilder, trailingBuilder, overrideProfilePicture = null,
}) {
  const history = useHistory();
  const navigateToMember = () => history.push(`/sequencer/users/${member._id}`);

  return (
    <InstitutionMemberListCard
      imageURL={overrideProfilePicture || member.avatarUrl}
      enforceImage
      title={`${member.firstName} ${member.lastName}`}
      description={member.email}
      nextToTitle={nextToNameBuilder ? nextToNameBuilder(member) : null}
      trailing={(trailingBuilder || (() => null))(member)}
      onClick={navigateToMember}
    />
  );
}

export default InstitutionMemberCard;
