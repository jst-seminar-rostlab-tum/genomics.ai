import React, { useState, useEffect } from 'react';
import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';
import PlusIcon from 'components/general/PlusIcon';
import TeamCard from 'components/teams/TeamCard';
import styles from './teamOverview.module.css';
import TeamService from 'shared/services/Team.service';
import TeamCreationDialog from 'components/teams/overview/TeamCreationDialog';

function TeamOverview() {
  const [teams, setTeams] = useState([]);
  useEffect(() => {
    TeamService.getMyTeams()
      .then(setTeams)
      .catch((e) => {
        console.error(e);
        if (e.response?.data) {
          alert(e.response.data);
        }
      });
  }, []);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  return (
    <>
      <Stack direction="row" className="stack" alignItems="Center">
        <Typography variant="h5" sx={{ pr: 1 }}>Your Teams</Typography>
        <PlusIcon onClick={() => setCreateOpen(true)} />
      </Stack>
      <br />
      {teams.length === 0 ? 'No teams.' : ''}
      {teams.map((team) => (
        <div key={team.id}>
          <TeamCard
            team={team}
            onLeft={(t) => onLeft(t)}
          />
          <div className={styles.cardSpacing} />
        </div>
      ))}
      <TeamCreationDialog
        open={createOpen}
        handleClose={() => setCreateOpen(false)}
        onCreated={(newTeam) => setTeams([...teams, newTeam])}
      />
    </>
  );
}

export default TeamOverview;
