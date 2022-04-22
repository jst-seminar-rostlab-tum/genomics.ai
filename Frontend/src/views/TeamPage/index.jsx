import React, { useState, useEffect } from 'react';
import {
  useParams,
} from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import JobList from 'components/teams/detail/TeamJobList';
import { getTeam } from 'shared/services/mock/teams';
import getUser from 'shared/services/mock/user';
import { getInstitution, queryIsAdminInstitutions } from 'shared/services/mock/institutions';
import Chip from '@mui/material/Chip';
import Stack from '@mui/material/Stack';
import TextField from '@mui/material/TextField';
import MenuItem from '@mui/material/MenuItem';
import Button from '@mui/material/Button';

export default function TeamPage({ sidebarShown }) {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user, setUser] = useState({});
  const [institution, setInstitution] = useState({});
  const [isMember, setIsMember] = useState(false);
  const [isAdmin, setIsAdmin] = useState(false);
  const [adminInstitutions, setAdminInstitutions] = useState([]);

  function updateIsMember() {
    setIsMember((team.memberIds || []).indexOf(user.id) > -1);
  }

  function updateIsAdmin() {
    setIsAdmin((team.adminIds || []).indexOf(user.id) > -1);
  }

  const handleDescriptionChange = (event) => {
    setTeam({
      ...team,
      description: event.target.value,
    });
  };

  useEffect(() => {
    getUser()
      .then((newUser) => { setUser(newUser); updateIsAdmin(); updateIsMember(); });
  }, [setUser, isAdmin, isMember]);

  useEffect(() => {
    getTeam(id)
      .then((newTeam) => { setTeam(newTeam); updateIsAdmin(); updateIsMember(); })
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam, isAdmin, isMember]);

  // Institution may be undefined
  useEffect(() => {
    getInstitution(team.institutionId)
      .then((newInstitution) => setInstitution(newInstitution));
  }, [team, setInstitution]);

  useEffect(() => {
    queryIsAdminInstitutions(user.id)
      .then((newAdminInstitutions) => setAdminInstitutions(newAdminInstitutions));
  }, [user, setAdminInstitutions]);

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title={team.name}
      rightOfTitle={(
        <HeaderOptions
          team={team}
          isAdmin={isAdmin}
          institution={institution}
          isMember={isMember}
          availableInstitutions={adminInstitutions}
          setInstitution={setInstitution}
        />
      )}
      replaceHeaderRight={
        (isAdmin && <AdminTeamHeaderRight team={team} setTeam={setTeam} />)
        || <UserTeamHeaderRight />
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
        <JobList teamId={id} forPart="geneMapper" />
      </section>
      <section>
        <h2>GeneCruncher</h2>
        <hr />
        <JobList teamId={id} forPart="geneCruncher" />
      </section>
      <section>
        <h2>Members</h2>
        <hr />
        TODO: implement
      </section>
    </HeaderView>
  );
}

function HeaderOptions({
  team, isAdmin, institution, availableInstitutions, setInstitution,
}) {
  const handleInstitutionChange = (event) => {
    setInstitution(event.target.value);
  };

  return (
    <Stack
      direction="row"
      alignItems="flex-end"
      spacing={2}
      sx={{ paddingLeft: '20px' }}
    >
      {(!isAdmin && institution) && <h4>{institution.name}</h4>}
      {!isAdmin && <Chip label={team.visibility} color="primary" />}

      {isAdmin
        && (
          <TextField
            id="select-institution"
            select
            label="Institution"
            value={institution}
            onChange={handleInstitutionChange}
            variant="standard"
          >
            {availableInstitutions.map((institutionOption) => (
              <MenuItem
                key={institutionOption.id}
                value={institutionOption}
              >
                {institutionOption.name}
              </MenuItem>
            ))}
          </TextField>
        )}

    </Stack>
  );
}

function AdminTeamHeaderRight({ team, setTeam }) {
  console.log(team);

  const updateVisibility = (event) => {
    setTeam({
      ...team,
      visibility: event.target.value,
    });
  };

  return (
    <TextField
      id="select-visibility"
      select
      label="Visibility"
      value={team.visibility}
      onChange={updateVisibility}
      variant="standard"
    >
      <MenuItem value="public">public</MenuItem>
      <MenuItem value="private">private</MenuItem>
      <MenuItem value="by institution">by institution</MenuItem>
    </TextField>
  );
}

function UserTeamHeaderRight() {
  return (
    <Button onClick={() => true}>Join</Button>
  );
}