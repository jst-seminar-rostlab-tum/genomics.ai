import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import styles from './teamMemberList.module.css';
import { useAuth } from 'shared/context/authContext';
import TeamService from 'shared/services/Team.service';

function TeamMemberList({ team, onMemberRemoved }) {
  const [user] = useAuth();
  const [members, setMembers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  useEffect(async () => {
    if (team.id == null) return;
    setIsLoading(true);
    setMembers(await TeamService.getMembers(team.id));
    setIsLoading(false);
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
        team.adminIds.includes(user.id) && user.id === member.id ? null : (
          <TeamMemberRemoveButton team={team} member={member} onRemoved={onMemberRemoved} />
        )
      )}
    />
  );
}

export default TeamMemberList;
