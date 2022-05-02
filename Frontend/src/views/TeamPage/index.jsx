import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import TeamJobList from 'components/teams/detail/TeamJobList';
import TeamMemberList from 'components/teams/detail/TeamMemberList';
import TeamAdminHeaderRight from 'components/teams/detail/TeamAdminHeaderRight';
import TeamUserHeaderRight from 'components/teams/detail/TeamUserHeaderRight';
import TeamHeaderOptions from 'components/teams/detail/TeamHeaderOptions';
import TeamInviteButton from 'components/teams/detail/TeamInviteButton';
import { getTeam } from 'shared/services/mock/teams';
import { getInstitution } from 'shared/services/mock/institutions';
import TextField from '@mui/material/TextField';
import { useAuth } from 'shared/context/authContext';

export default function TeamPage({ sidebarShown }) {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user] = useAuth();
  const [institution, setInstitution] = useState({});

  const handleDescriptionChange = (event) => {
    setTeam({
      ...team,
      description: event.target.value,
    });
  };

  useEffect(() => {
    getTeam(parseInt(id, 10))
      .then(setTeam)
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam]);

  // Institution may be undefined
  useEffect(() => {
    getInstitution(team.institutionId)
      .then((newInstitution) => setInstitution(newInstitution));
  }, [team, setInstitution]);

  const isAdmin = team.adminIds ? team.adminIds.includes(user.id) : false;

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title={team.name}
      rightOfTitle={(
        <TeamHeaderOptions
          team={team}
          isAdmin={isAdmin}
          institution={institution}
        />
      )}
      replaceHeaderRight={
        (isAdmin && <TeamAdminHeaderRight team={team} setTeam={setTeam} />)
        || <TeamUserHeaderRight institution={institution} team={team} user={user} />
      }
    >
      <br />
      <section>
        <h2>Description</h2>
        <hr />
        <TextField
          id="description"
          multiline
          minRows={3}
          maxRows={5}
          value={team.description}
          InputProps={{
            readOnly: !isAdmin,
          }}
          fullWidth
          onChange={handleDescriptionChange}
          variant="standard"
        />
      </section>
      <section>
        <h2>GeneMapper</h2>
        <hr />
        <TeamJobList team={team} forPart="geneMapper" />
      </section>
      <section>
        <h2>GeneCruncher</h2>
        <hr />
        <TeamJobList team={team} forPart="geneCruncher" />
      </section>
      <section>
        <h2>Members</h2>
        <hr />
        <TeamMemberList team={team} />
      </section>
      {isAdmin && <TeamInviteButton team={team} />}
    </HeaderView>
  );
}
