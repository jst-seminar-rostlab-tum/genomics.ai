import React, { useState, useEffect } from 'react';
import Button from '@mui/material/Button';
import TeamLeaveButton from 'components/teams/overview/TeamLeaveButton';

function TeamUserHeaderRight({ team, user }) {
  const [isMember, setIsMember] = useState(false);

  function updateIsMember() {
    setIsMember((team.memberIds || []).includes(user.id));
  }

  const onLeft = () => {
    setIsMember(false);
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
    <Button onClick={() => true}>Join</Button>
  );
}

export default TeamUserHeaderRight;
