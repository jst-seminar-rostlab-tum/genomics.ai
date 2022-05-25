import React, { useState, useEffect } from 'react';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';
import TeamJoinButton from 'components/teams/detail/TeamJoinButton';
import { Tooltip } from '@mui/material';

function TeamUserHeaderRight({
  institution, team, user, updateTeam,
}) {
  const [isMember, setIsMember] = useState(false);

  function updateIsMember() {
    setIsMember((team.memberIds || []).includes(user._id));
  }

  const onLeft = () => {
    updateTeam();
  };

  const onJoin = () => {
    updateTeam();
  };

  useEffect(() => {
    updateIsMember();
  }, [team]);

  if (isMember) {
    return (
      <TeamLeaveButton team={team} onLeft={onLeft} />
    );
  }

  const canJoin = (team.visibility === 'PUBLIC') || (team.visibility === 'BY_INSTITUTION' && institution.memberIds.includes(user._id))
    || (team.invitedMemberIds.includes(user._id));

  return (
    <Tooltip title={canJoin ? '' : `This team is ${team.visibility.toLowerCase()} and you haven't been invited${team.visibility === 'BY_INSTITUTION' ? " or you're a member of this team's institution" : ''}.`}>
      <div>
        <TeamJoinButton
          isDisabled={!canJoin}
          team={team}
          onJoin={onJoin}
        />
      </div>
    </Tooltip>
  );
}

export default TeamUserHeaderRight;
