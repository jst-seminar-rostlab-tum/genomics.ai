import React, { useState, useEffect } from 'react';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import TeamJoinButton from 'components/teams/detail/TeamJoinButton';

function TeamUserHeaderRight({ institution, team, user }) {
  const [isMember, setIsMember] = useState(false);

  function updateIsMember() {
    setIsMember((team.memberIds || []).includes(user.id));
  }

  const onLeft = () => {
    setIsMember(false);
  };

  const onJoin = () => {
    setIsMember(true);
  };

  useEffect(() => {
    updateIsMember();
  }, [team]);

  if (isMember) {
    return (
      <TeamLeaveButton team={team} onLeft={onLeft} />
    );
  }

  return (
    <TeamJoinButton
      isDisabled={(team.visibility === 'private' && (!team.invitedMemberIds.includes(user.id)))
      || (team.visibility === 'by institution' && (!institution.memberIds.includes(user.id) && !institution.adminIds.includes(user.id)))}
      team={team}
      onJoin={onJoin}
    />
  );
}

export default TeamUserHeaderRight;
