import React, { useState, useEffect } from 'react';
import { CircularProgress, Stack } from '@mui/material';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import TeamMemberMakeAdminButton from '../TeamMemberMakeAdminButton';
import styles from './teamMemberList.module.css';
import { useAuth } from 'shared/context/authContext';
import TeamService from 'shared/services/Team.service';

function TeamMemberList({
  team, updateTeam,
}) {
  const [user] = useAuth();
  const [members, setMembers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  useEffect(() => {
    if (team.id == null) return;
    setIsLoading(true);
    TeamService.getMembers(team.id)
      .then((newMembers) => {
        setMembers(newMembers);
        setIsLoading(false);
      });
  }, [team]);

  if (isLoading) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      members={members}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {team.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        team.adminIds.includes(user._id) && user._id !== member.id ? (
          <Stack direction="row" spacing={1}>
            <TeamMemberMakeAdminButton
              team={team}
              member={member}
              updateTeam={updateTeam}
            />
            <TeamMemberRemoveButton
              team={team}
              member={member}
              updateTeam={updateTeam}
            />
          </Stack>
        ) : null
      )}
    />
  );
}

export default TeamMemberList;
