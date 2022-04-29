import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import getProfile from 'shared/services/profile';
import styles from './teamMemberList.module.css';

function TeamMemberList({ team }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getProfile().then(setUser);
  }, [setUser]);

  if (!team.adminIds?.length || !team.memberIds?.length) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      memberIds={[...team.adminIds, ...team.memberIds]}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {team.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        team.adminIds.includes(user.id) && user.id === member.id ? null : (
          <TeamMemberRemoveButton team={team} member={member} />
        )
      )}
    />
  );
}

export default TeamMemberList;
